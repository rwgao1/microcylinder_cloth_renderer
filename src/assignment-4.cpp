#include "material.h"
#include "muni/camera.h"
#include "muni/common.h"
#include "muni/image.h"
#include "muni/material.h"
#include "muni/math_helpers.h"
#include "muni/obj_loader.h"
#include "muni/ray_tracer.h"
#include "muni/sampler.h"
#include "muni/scenes/box.h"
#include "muni/triangle.h"
#include "ray_tracer.h"
#include "spdlog/spdlog.h"
#include "triangle.h"
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace muni;

RayTracer::Octree octree{};

std::tuple<Vec3f, float> sampleHemi(Vec3f normal, Vec2f u) {
    // =============================================================================================
    float theta = std::acos(1.0f - 2.0f * u[0]);
    float phi = 2 * M_PI * u[1];

    Vec3f dir = Vec3f{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};

    if (dot(dir, normal) < 0) {
        dir = -dir;
    }

    return {normalize(dir), 1.0f / (2.0f * M_PI)};
    // =============================================================================================
}

std::tuple<Vec3f, float> sampleHemiCos(Vec3f normal, Vec2f u) {
    // =============================================================================================
    float theta = std::acos(sqrt(u[0]));
    float sin_theta = std::sin(theta);
    float cos_theta = sqrt(u[0]);

    float phi = 2 * M_PI * u[1];
    float sin_phi = std::sin(phi);
    float cos_phi = std::cos(phi);

    Vec3f dir = Vec3f{sin_theta * cos_phi, sin_theta * sin_phi, cos_theta};

    dir = from_local(dir, normal);

    return {normalize(dir), cos(theta) / M_PI};
    // =============================================================================================
}

/** Offset the ray origin to avoid self-intersection.
    \param[in] ray_pos The original ray origin.
    \param[in] normal The normal of the surface at the hit point.
    \return The offset ray origin.
*/
Vec3f offset_ray_origin(Vec3f ray_pos, Vec3f normal) {
    return ray_pos + EPS * normal;
}

/** Check if the triangle is an emitter.
    \param[in] tri The triangle to check
    \return True if the triangle is an emitter, false otherwise.
*/
bool is_emitter(const Triangle &tri) { return tri.emission != Vec3f{0.0f}; }

/** Evaluate the radiance of the area light. We **do not** check whether the hit
 point is on the light source, so make sure
 *  the hit point is on the light source before calling this function.
    \param[in] light_dir The **outgoing** direction from the light source to the
 scene. \return The radiance of the light source.
*/
Vec3f eval_area_light(const Vec3f light_dir) {
    if (dot(light_dir, BoxScene::light_normal) > 0.0f)
        return BoxScene::light_color;
    return Vec3f{0.0f};
}

/** Sample a point on the area light with a uniform distribution.
    \param[in] samples A 2D uniform random sample.
    \return A tuple containing the sampled position, the normal of the light
 source, and the PDF value.
*/
std::tuple<Vec3f, Vec3f, float> sample_area_light(Vec2f samples) {
    // =============================================================================================
    // TODO: Implement this function
    // =============================================================================================
    Vec3f light_pos{BoxScene::light_x, BoxScene::light_y, BoxScene::light_z};
    Vec3f light_right{BoxScene::light_len_x, 0.0f, 0.0f};
    Vec3f light_up{0.0f, BoxScene::light_len_y, 0.0f};

    Vec3f xp = light_pos + (samples[0] * light_right) + (samples[1] * light_up);

    return {xp, BoxScene::light_normal, BoxScene::inv_light_area};
}


Vec3f shade_with_light_sampling(Triangle tri, Vec3f p, Vec3f wo) {
    // =============================================================================================
    // TODO: Implement this function
    // Please refer to lecture 9, page 20&21 for the details of the implementation.
    std::tuple<Vec3f, float> sample = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]) ? 
        std::get<Lambertian>(BoxScene::materials[tri.material_id]).sample(tri.face_normal, UniformSampler::next2d()) :
        std::get<Microfacet>(BoxScene::materials[tri.material_id]).sample(wo, tri.face_normal, UniformSampler::next2d());

    auto [wi, pdf] = sample;

    Vec3f f_r;
    if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])) {
        f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
    } else if (std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id])) {
        f_r = std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, wi, tri.face_normal);

    } else {
        f_r = Vec3f{0.0f};
    }

    // Contribution from the light source
    Vec3f L_dir{0.0f};
    // Uniformly sample the light at x
    auto [xp, light_normal, pdf_light] = sample_area_light(UniformSampler::next2d());
    // Shoot a ray from p to x points from light to point
    Vec3f wi_light = -normalize(xp - p);
    
    Vec3f f_r_light = std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]) ? 
        std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, -wi_light, tri.face_normal) : // check direction
        f_r;

    pdf_light = length_squared(xp - p) * BoxScene::inv_light_area / std::max<float>(EPS, dot(light_normal, wi_light));

    auto [light_is_ray_hit, light_t_min, light_nearest_tri] = RayTracer::closest_hit(p, -wi_light, octree, BoxScene::triangles);
    // If the ray is not blocked in the middle
    if (light_is_ray_hit && is_emitter(light_nearest_tri)) {
        Vec3f L_light = eval_area_light(wi_light);
        // Compute the contribution from the light source
        L_dir = f_r_light * L_light * std::max<float>(0.0f, dot(-wi_light, tri.face_normal)) / std::max(EPS, pdf_light);
    }
    // Contribution from other reflectors
    Vec3f L_indir{0.0f};

    // Test Russian Roulette with probability p_rr = 0.8f
    const float p_rr = 0.8f;

    if (UniformSampler::next1d() > p_rr) {
        return L_dir;
    }
    

    // Trace the new ray
    auto [is_ray_hit, t_min, nearest_tri] = RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);

    // If the ray hit a non-emitting object at q
    if (!is_emitter(nearest_tri)) {
        // Vec3f q = offset_ray_origin(p + t_min * wi, nearest_tri.face_normal);
        Vec3f q = offset_ray_origin(p + t_min * wi, nearest_tri.face_normal);
        L_indir = f_r * shade_with_light_sampling(nearest_tri, q, -wi) * std::max(0.0f, dot(wi, tri.face_normal)) / pdf / p_rr;
    }

    return L_dir + L_indir;
    // =============================================================================================
}

Vec3f shade_cloth(Triangle tri, Vec3f p, Vec3f wo);
Vec3f shade_thread(Triangle tri, Vec3f p, Vec3f wr, Vec3f t, Vec3f wi, float pdf, std::vector<Vec3f> ts_offset, std::vector<Vec3f> ns_offset, float W, int direction, int N1, int N2, float a1, float a2);

Vec3f shade_thread(Triangle tri, Vec3f p, Vec3f wr, Vec3f t, Vec3f wi, float pdf_wi, std::vector<Vec3f> ts_offset, std::vector<Vec3f> ns_offset, float W, int direction, int N1, int N2, float a1, float a2) {
    Vec3f surface_normal = tri.face_normal;
    float Nj = ts_offset.size();
    Vec3f L_rj{0.0f};
    float P_total_dir = 0.0f;
    float P_total_indir = 0.0f;
    float P_total = 0.0f;

    for (int i = 0; i < Nj; i++) {
        // direction 0 is x thread, 1 is y
        Vec3f t_offset = ts_offset.at(i);
        Vec3f local_normal_offset = ns_offset.at(i);

        // auto [wi, pdf_wi] = std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id]) ? 
        //     std::get<Lambertian>(BoxScene::materials[tri.material_id]).sample(local_normal_offset, UniformSampler::next2d()) :
        //     std::get<Cloth>(BoxScene::materials[tri.material_id]).sample(wr, local_normal_offset, t_offset, UniformSampler::next3d());

       

        Vec3f L_dir{0.0f};
        Vec3f L_indir{0.0f};

        // Uniformly sample the light at x
        auto [xp, light_normal, pdf_light] = sample_area_light(UniformSampler::next2d());
        // Shoot a ray from p to x points from light to point
        Vec3f wi_light = -normalize(xp - p);

        Vec3f f_r_light;
        Vec3f f_r;

        float P_dir = 1.0f;
        float P_indir = 1.0f;


        pdf_light = length_squared(xp - p) * BoxScene::inv_light_area / safe(dot(light_normal, wi_light));

        auto [light_is_ray_hit, light_t_min, light_nearest_tri] = RayTracer::closest_hit(p, -wi_light, octree, BoxScene::triangles);
        // If the ray is not blocked in the middle
        if (light_is_ray_hit && is_emitter(light_nearest_tri)) {

            if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])) {
                f_r_light = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
            } else if (std::holds_alternative<Cloth>(BoxScene::materials[tri.material_id])) {
                f_r_light = std::get<Cloth>(BoxScene::materials[tri.material_id]).eval(-wi_light, wr, local_normal_offset, t_offset, direction);
                // P_dir = std::get<Cloth>(BoxScene::materials[tri.material_id]).P(-wi_light, wr, t_offset, local_normal_offset);
            } else {
                f_r_light = Vec3f{0.0f};
            }

            Vec3f L_light = eval_area_light(wi_light);
            // Compute the contribution from the light source
            L_dir = P_dir * f_r_light * L_light * std::max<float>(0.0f, dot(-wi_light, local_normal_offset)) / std::max(EPS, pdf_light);
        }

        L_rj += L_dir;

        // Test RR
        const float p_rr = 0.15f;
        if (UniformSampler::next1d() > p_rr) {
            P_total += P_dir;
            continue;
        }

        // Trace ray
        auto [is_ray_hit, t_min, nearest_tri] = RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);

        if (is_ray_hit && !is_emitter(nearest_tri)) {

            if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])) {
                f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
            } else if (std::holds_alternative<Cloth>(BoxScene::materials[tri.material_id])) {
                f_r = std::get<Cloth>(BoxScene::materials[tri.material_id]).eval(wi, wr, local_normal_offset, t_offset, direction);
                // P_indir = std::get<Cloth>(BoxScene::materials[tri.material_id]).P(wi, wr, t_offset, local_normal_offset);
            } else {
                f_r = Vec3f{0.0f};
            }
            // TODO: offset ray origin might be jank
            Vec3f q = offset_ray_origin(p + t_min * wi, nearest_tri.face_normal);
            L_indir += P_indir * f_r * shade_cloth(nearest_tri, q, -wi) * std::max(0.0f, dot(wi, local_normal_offset)) / safe(pdf_wi) / safe(p_rr);
        }
        L_rj += L_indir;

        P_total += P_dir + P_indir;
    }

    return L_rj;
    // float Q = (a1 / N1) * P_total + (a2 / N2) * P_total + (1 - a1 - a2) * std::max<float>(0.0f, dot(wr, surface_normal));
    return L_rj / Nj;
}

Vec3f shade_cloth(Triangle tri, Vec3f p, Vec3f wo) {

    if (std::holds_alternative<Cloth>(BoxScene::materials[tri.material_id])) {
        int N1 = std::get<Cloth>(BoxScene::materials[tri.material_id]).N1;
        int N2 = std::get<Cloth>(BoxScene::materials[tri.material_id]).N2;
        Vec3f t1 = std::get<Cloth>(BoxScene::materials[tri.material_id]).t1;
        Vec3f t2 = std::get<Cloth>(BoxScene::materials[tri.material_id]).t2;
        float a1 = std::get<Cloth>(BoxScene::materials[tri.material_id]).a1;
        float a2 = std::get<Cloth>(BoxScene::materials[tri.material_id]).a2;

        Vec3f surface_n = tri.face_normal;
        Vec3f surface_n_xz = normalize(Vec3f{surface_n.x, 0.0f, surface_n.z});
        Vec3f surface_n_yz = normalize(Vec3f{0.0f, surface_n.y, surface_n.z});
        
        // get orthogonal t's
        t1 = normalize(t1 - project(t1, surface_n_xz));
        t2 = normalize(t2 - project(t2, surface_n_yz));

        auto [alpha1s, alpha2s] = std::get<Cloth>(BoxScene::materials[tri.material_id]).sample_offsets(N1, N2);
        auto [wi, pdf_wi] = sampleHemiCos(tri.face_normal , UniformSampler::next2d());

        float W = 0.0f;
        std::vector<Vec3f> t1s_offset;
        std::vector<Vec3f> t2s_offset;
        std::vector<Vec3f> n1s_offset;
        std::vector<Vec3f> n2s_offset;
        for (float alpha1 : alpha1s) {
            t1s_offset.push_back(rotate_y(t1, alpha1));
            n1s_offset.push_back(rotate_y(surface_n, alpha1));
            W += std::get<Cloth>(BoxScene::materials[tri.material_id]).P(wi, wo, t1s_offset.back(), n1s_offset.back());
        }

        for (float alpha2 : alpha2s) {
            t2s_offset.push_back(rotate_x(t2, alpha2));
            n2s_offset.push_back(rotate_x(surface_n, alpha2));
            W += std::get<Cloth>(BoxScene::materials[tri.material_id]).P(wi, wo, t2s_offset.back(), n2s_offset.back());
        }

        Vec3f L_r1 = shade_thread(tri, p, wo, t1, wi, pdf_wi, t1s_offset, n1s_offset, W, 0, N1, N2, a1, a2);
        Vec3f L_r2 = shade_thread(tri, p, wo, t2, wi, pdf_wi, t2s_offset, n2s_offset, W, 1, N1, N2, a1, a2);

        return a1 * L_r1 + a2 * L_r2;
    } else {
        // Lambertian

        std::tuple<Vec3f, float> sample = sampleHemi(tri.face_normal, UniformSampler::next2d());

        auto [wi, pdf] = sample;

        Vec3f f_r;
        if (std::holds_alternative<Lambertian>(BoxScene::materials[tri.material_id])) {
            f_r = std::get<Lambertian>(BoxScene::materials[tri.material_id]).eval();
        } else {
            f_r = Vec3f{0.0f};
        }

        // Contribution from the light source
        Vec3f L_dir{0.0f};
        // Uniformly sample the light at x
        auto [xp, light_normal, pdf_light] = sample_area_light(UniformSampler::next2d());
        // Shoot a ray from p to x points from light to point
        Vec3f wi_light = -normalize(xp - p);
        
        Vec3f f_r_light = std::holds_alternative<Microfacet>(BoxScene::materials[tri.material_id]) ? 
            std::get<Microfacet>(BoxScene::materials[tri.material_id]).eval(wo, -wi_light, tri.face_normal) : // check direction
            f_r;

        pdf_light = length_squared(xp - p) * BoxScene::inv_light_area / std::max<float>(EPS, dot(light_normal, wi_light));

        auto [light_is_ray_hit, light_t_min, light_nearest_tri] = RayTracer::closest_hit(p, -wi_light, octree, BoxScene::triangles);
        // If the ray is not blocked in the middle
        if (light_is_ray_hit && is_emitter(light_nearest_tri)) {
            Vec3f L_light = eval_area_light(wi_light);
            // Compute the contribution from the light source
            L_dir = f_r_light * L_light * std::max<float>(0.0f, dot(-wi_light, tri.face_normal)) / std::max(EPS, pdf_light);
        }
        // Contribution from other reflectors
        Vec3f L_indir{0.0f};

        // Test Russian Roulette with probability p_rr = 0.8f
        const float p_rr = 0.15f;

        if (UniformSampler::next1d() > p_rr) {
            return L_dir;
        }
        

        // Trace the new ray
        auto [is_ray_hit, t_min, nearest_tri] = RayTracer::closest_hit(p, wi, octree, BoxScene::triangles);

        // If the ray hit a non-emitting object at q
        if (!is_emitter(nearest_tri)) {
            // Vec3f q = offset_ray_origin(p + t_min * wi, nearest_tri.face_normal);
            Vec3f q = offset_ray_origin(p + t_min * wi, nearest_tri.face_normal);
            L_indir = f_r * shade_cloth(nearest_tri, q, -wi) * std::max(0.0f, dot(wi, tri.face_normal)) / safe(pdf) / safe(p_rr);
        }

        return L_dir + L_indir;
    }
}


Vec3f path_tracing_with_light_sampling(Vec3f ray_pos, Vec3f ray_dir) {
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit) return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri)) return eval_area_light(-ray_dir);

    return shade_cloth(nearest_tri, hit_position, -ray_dir);
}


void test_conversions() {
    Vec3f v{-2.0f, 7.0f, 4.0f};
    auto [phi, theta] = to_azimuthal_longitudinal(v, Vec3f{1.0f, 0.0f, 0.0f});
    spdlog::info("phi: {}, theta: {}", phi, theta);
    Vec3f vp = to_cartesian(phi, theta, Vec3f{1.0f, 0.0f, 0.0f});

    spdlog::info("vp: {}", vp);
    spdlog::info("v: {}", normalize(v));    
}


void test_get_default_t() {
    Vec3f t1{1.0f, 0.0f, 0.0f};
    Vec3f t2{0.0f, 1.0f, 0.0f};

    Vec3f surface_n{-1.0f, 1.0f, 1.0f};
    Vec3f surface_n_xz = normalize(Vec3f{surface_n.x, 0.0f, surface_n.z});
    Vec3f surface_n_yz = normalize(Vec3f{0.0f, surface_n.y, surface_n.z});

    t1 = normalize(t1 - project(t1, surface_n_xz));
    t2 = normalize(t2 - project(t2, surface_n_yz));

    spdlog::info("t1: {}, t2: {}", t1, t2);

}

void add_to_scene(std::string obj_path, int material_id) {
    std::vector<Triangle> obj_triangles = load_obj(obj_path, material_id);
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
}


int main(int argc, char **argv) {

    spdlog::info("\n"
                 "----------------------------------------------\n"
                 "Welcome to CS 190I Assignment 4: Microfacet Materials\n"
                 "----------------------------------------------");
    const unsigned int max_spp = 1024;
    const unsigned int image_width = 1080;
    const unsigned int image_height = 1080;
    // Some prepereations
    Image image{.width = image_width,
                .height = image_height,
                .pixels = std::vector<Vec3f>(image_width * image_height)};
    Camera camera{.vertical_field_of_view = 38.6f,
                  .aspect = static_cast<float>(image_width) / image_height,
                  .focal_distance = 0.8f,
                  .position = Vec3f{0.278f, 0.8f, 0.2744f},
                  .view_direction = Vec3f{0.0f, -1.0f, 0.0f},
                  .up_direction = Vec3f{0.0f, 0.0f, 1.0f},
                  .right_direction = Vec3f{-1.0f, 0.0f, 0.0f}};
    camera.init(); 
    UniformSampler::init(190);


    // =============================================================================================


    add_to_scene("./draped.obj", 9);
    add_to_scene("./right.obj", 9);


    octree.build_octree(BoxScene::triangles);
    // =============================================================================================
    // Path Tracing with light sampling
    spdlog::info("Path Tracing with light sampling: rendering started!");

    omp_set_num_threads(12);
    #pragma omp parallel for

    for (int y = 0; y < image.height; y++) {
        if (y % 50 == 0) {
            spdlog::info("Rendering row {} / {} \r", y, image.height);
        }
        for (int x = 0; x < image.width; x++) {
            image(x, y) = Vec3f{0.0f};
            for (int sample = 0; sample < max_spp; sample++) {
                const float u = (x + UniformSampler::next1d()) / image.width;
                const float v = (y + UniformSampler::next1d()) / image.height;
                Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
                image(x, y) += clamp(path_tracing_with_light_sampling(
                                         camera.position, ray_direction),
                                     Vec3f(0.0f), Vec3f(50.0f));
            }
            image(x, y) /= (float)max_spp;
        }
    }
    spdlog::info("Path Tracing with light sampling: Rendering finished!");
    // image.save_with_tonemapping("./cloth" + std::to_string(left_pillow_material_id) + "_" + std::to_string(right_pillow_material_id) + ".png");
    image.save_with_tonemapping("./scene.png");


    // =============================================================================================
    // Path Tracing with MIS
    // spdlog::info("Path Tracing with MIS: Rendering started!");
    // for (int y = 0; y < image.height; y++) {
    //     if (y % 50 == 0) {
    //         spdlog::info("Rendering row {} / {} \r", y, image.height);
    //     }
    //     for (int x = 0; x < image.width; x++) {
    //         image(x, y) = Vec3f{0.0f};
    //         for (int sample = 0; sample < max_spp; sample++) {
    //             const float u = (x + UniformSampler::next1d()) / image.width;
    //             const float v = (y + UniformSampler::next1d()) / image.height;
    //             Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
    //             image(x, y) +=
    //                 clamp(path_tracing_with_MIS(camera.position, ray_direction),
    //                       Vec3f(0.0f), Vec3f(50.0f));
    //         }
    //         image(x, y) /= (float)max_spp;
    //     }
    // }
    // spdlog::info("Path Tracing with MIS: Rendering finished!");
    // image.save_with_tonemapping("./path_tracing_with_MIS.png");

    // =============================================================================================
    return 0;
}

