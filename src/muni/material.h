#pragma once
#include "common.h"
#include "math_helpers.h"
#include "sampler.h"
#include <algorithm>
#include <cmath>


namespace muni {

struct Lambertian {
    Vec3f albedo;

    /** Evaluates the BRDF for the Lambertian material.
      \return The BRDF (fr) value.
  */
    Vec3f eval() const { 
        // =============================================================================================
        // TODO: Implement this function
        return albedo / M_PI;
        // =============================================================================================
     }

    /** Samples the BRDF for the Lambertian material.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
        // =============================================================================================
        // TODO: Implement this function
        float theta = std::acos(1.0f - 2.0f * u[0]);
        float phi = 2 * M_PI * u[1];

        Vec3f dir = Vec3f{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};

        if (dot(dir, normal) < 0) {
            dir = -dir;
        }

        return {normalize(dir), 1.0f / (2.0f * M_PI)};
        // =============================================================================================
    }
    /** Computes the PDF for the Lambertian material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        // =============================================================================================
        // TODO: Implement this function
        return 1.0f / (2.0f * M_PI);
        // =============================================================================================
    }

};

struct Cloth {
  Vec3f A1;
  Vec3f A2;
  float gamma_s1;
  float gamma_s2;
  float gamma_v1;
  float gamma_v2;
  float k_d1;
  float k_d2;
  Vec3f eta3;
  float a1;
  float a2;
  Vec3f t1;
  Vec3f t2;
  std::vector<float> t1_offsets;
  std::vector<float> t1_offset_lengths;
  std::vector<float> t2_offsets;
  std::vector<float> t2_offset_lengths;
  int N1;
  int N2;

  std::tuple<std::vector<float>, std::vector<float> > sample_offsets(int N1, int N2) const {
    std::vector<float> offsets1;
    std::vector<float> offsets2;

    for (int i = 0; i < N1; i++) {
      offsets1.push_back(sample_t1_offset(UniformSampler::next1d()));
    }

    for (int i = 0; i < N2; i++) {
      offsets2.push_back(sample_t2_offset(UniformSampler::next1d()));
    }

    return {offsets1, offsets2};
  }

  float sample_t1_offset(float u) const {
    // inverse sample offset_lengths
    float sum = 0;
    for (int i = 0; i < t1_offset_lengths.size(); i++) {
      sum += t1_offset_lengths.at(i);
    }

    std::vector<float> cdf = t1_offset_lengths;
    for (int i = 1; i < cdf.size(); i++) {
      cdf.at(i) += cdf.at(i - 1);
    }

    for (int i = 0; i < cdf.size(); i++) {
      cdf.at(i) /= sum;
    }

    auto low = std::lower_bound(cdf.begin(), cdf.end(), u);
    int index = std::distance(cdf.begin(), low);

    float percent_between;
    if (t1_offsets.size() == 1) {
      percent_between = u;
    } else {
      if (index == 0) {
        percent_between = u;
      } else {
        percent_between = (u - cdf.at(index - 1)) / (cdf.at(index) - cdf.at(index - 1));
      }
    }
    
    return M_PI / 180 * (t1_offsets.at(index) + percent_between * (t1_offsets.at(index + 1) - t1_offsets.at(index)));
  }

  float sample_t2_offset(float u) const {
    // inverse sample offset_lengths
    float sum = 0;
    for (int i = 0; i < t2_offset_lengths.size(); i++) {
      sum += t2_offset_lengths.at(i);
    }

    std::vector<float> cdf = t2_offset_lengths;
    for (int i = 1; i < cdf.size(); i++) {
      cdf.at(i) += cdf.at(i - 1);
    }

    for (int i = 0; i < cdf.size(); i++) {
      cdf.at(i) /= sum;
    }

    auto low = std::lower_bound(cdf.begin(), cdf.end(), u);
    int index = std::distance(cdf.begin(), low);

    float percent_between;
    if (t2_offsets.size() == 1) {
      percent_between = u;
    } else {
      if (index == 0) {
        percent_between = u;
      } else {
        percent_between = (u - cdf.at(index - 1)) / (cdf.at(index) - cdf.at(index - 1));
      }
    }
    
    return M_PI / 180 * (t2_offsets.at(index) + percent_between * (t2_offsets.at(index + 1) - t2_offsets.at(index)));
  }


  float M(Vec3f w_i, Vec3f w_r, Vec3f t, Vec3f normal) const {
    Vec3f wi_local = to_local(w_i, normal);
    Vec3f wr_local = to_local(w_r, normal);
    Vec3f t_local = to_local(t, normal);

    auto [phi_i, theta_i] = to_azimuthal_longitudinal(wi_local, t_local);
    auto [phi_r, theta_r] = to_azimuthal_longitudinal(wr_local, t_local);

    float phi_d = phi_i - phi_r;
    float sigma = 20 * M_PI / 180; // between 15-25 deg
    float u = gaussian_u(phi_d, 0, sigma);
    float M_wi = std::max(0.0f, cos(phi_i));
    float M_wr = std::max(0.0f, cos(phi_r));

    return (1 - u) * M_wi * M_wr + u * std::min(M_wi, M_wr);
  }


  float P(Vec3f w_i, Vec3f w_r, Vec3f t, Vec3f normal) const {
    Vec3f wi_local = to_local(w_i, normal);
    Vec3f wr_local = to_local(w_r, normal);
    Vec3f t_local = to_local(t, normal);

    float psi_i = get_psi(wi_local, t_local);
    float psi_r = get_psi(wr_local, t_local);
    float psi_d = psi_i - psi_r;

    float sigma = 20 * M_PI / 180; // between 15-25 deg
    float u = gaussian_u(psi_d, 0, sigma);

    float P_wi = std::max(0.0f, cos(psi_i));
    float P_wr = std::max(0.0f, cos(psi_r));

    return (1 - u) * P_wi * P_wr + u * std::min(P_wi, P_wr);
  }

  float gaussian(float x, float mu, float sigma) const {
    return exp(-pow(x - mu, 2) / pow(sigma, 2));
  }

  float gaussian_u(float x, float mu, float sigma) const {
    return gaussian(x, mu, sigma) / (sigma * sqrt(M_PI));
  }

  Vec3f F_r_dielectric(float cos_theta) const {
    Vec3f temp = (eta3 - 1) / (eta3 + 1);
    Vec3f R0{temp[0] * temp[0], temp[1] * temp[1], temp[2] * temp[2]};

    return R0 + (1 - R0) * (float) pow(1 - cos_theta, 5);
  }

  Vec3f F_t_dielectric(float theta) const {
    return 1 - F_r_dielectric(theta);
  }

  // surface scattering
  Vec3f f_rs(Vec3f wi, Vec3f wr, Vec3f t, int dir) const {
    float gamma_s = dir == 0 ? gamma_s1 : gamma_s2;
    auto [phi_i, theta_i] = to_azimuthal_longitudinal(wi, t);
    auto [phi_r, theta_r] = to_azimuthal_longitudinal(wr, t);

    float phi_d = phi_i - phi_r;
    float theta_d = (theta_i - theta_r) / 2;
    float theta_h = (theta_i + theta_r) / 2;

    float g = gaussian_u(theta_h, 0.0f, gamma_s);

    Vec3f F = F_r_dielectric(cos(theta_d) * cos(phi_d / 2));

    return F * cos(phi_d / 2) * g;
  }

  // volume scattering
  Vec3f f_rv(Vec3f wi, Vec3f wr, Vec3f t, int dir) const {
    // assume local space
    float gamma_v = dir == 0 ? gamma_v1 : gamma_v2;
    float k_d = dir == 0 ? k_d1 : k_d2;
    Vec3f A = dir == 0 ? A1 : A2;

    auto [phi_i, theta_i] = to_azimuthal_longitudinal(wi, t);
    auto [phi_r, theta_r] = to_azimuthal_longitudinal(wr, t);

    float theta_h = (theta_i + theta_r) / 2;

    Vec3f F_ti = 1 - F_t_dielectric(wi.z);
    Vec3f F_tr = 1 - F_t_dielectric(wr.z);

    Vec3f F = F_ti * F_tr;

    float g = gaussian_u(theta_h, 0, gamma_v);

    return A * (F * (1 - k_d) * g + k_d) / (cos(theta_i) * cos(theta_r));
  }
  
  Vec3f eval(Vec3f wi, Vec3f wr, Vec3f normal, Vec3f t, int dir) const {
    Vec3f wi_local = to_local(wi, normal);
    Vec3f wr_local = to_local(wr, normal);
    Vec3f t_local = to_local(t, normal);

    auto [phi_i, theta_i] = to_azimuthal_longitudinal(wi_local, t_local);
    auto [phi_r, theta_r] = to_azimuthal_longitudinal(wr_local, t_local);

    float theta_d = (theta_i - theta_r) / 2;

    Vec3f frs = f_rs(wi_local, wr_local, t_local, dir);
    Vec3f frv = f_rv(wi_local, wr_local, t_local, dir);

    return M(wi, wr, t, normal) * (frs + frv) / (float) pow(cos(theta_d), 2);
  }

  std::tuple<Vec3f, float> sample(Vec3f wr_world, Vec3f normal, Vec3f t, Vec3f u, int dir) const {
    float gamma_s = dir == 0 ? gamma_s1 : gamma_s2;
    float gamma_v = dir == 0 ? gamma_v1 : gamma_v2;
    float k_d = dir == 0 ? k_d1 : k_d2;

    Vec3f wo_local = to_local(wr_world, normal);
    Vec3f t_local = to_local(t, normal);

    auto [phi_r, theta_r] = to_azimuthal_longitudinal(wo_local, t_local);
    // with probability 1-k_d sample surface scattering
    // with probability k_d sample volume scattering
    float theta_i;
    float phi_i;

    float A = atan((theta_r + M_PI_2) / (2 * gamma_s));
    float B = atan((theta_r - M_PI_2) / (2 * gamma_s));

    if (std::isinf(A) || std::isinf(B)) {
      spdlog::info("A: {}", A);
      spdlog::info("B: {}", B);
      exit(0);
    }
    

    if (u[0] < k_d) {
      theta_i = acos(2 * u[1] - 1);

    } else {
      theta_i = 2 * gamma_s * tan(u[1] * (A - B) + B) - theta_r;
    }

    phi_i = phi_r - 2 * asin(2 * u[2] - 1);

    float theta_h = (theta_i + theta_r) / 2;
    float phi_d = phi_i - phi_r;

    float pdf_theta_i_v = cos(theta_i) / 2;
    float pdf_theta_i_s = 1.0f / (2 * (A - B)) * gamma_s / safe(theta_h * theta_h + gamma_s * gamma_s);

    float pdf_theta_i = k_d * pdf_theta_i_v + (1 - k_d) * pdf_theta_i_s;
    float pdf_phi_i = 0.25 * cos(phi_d / 2);

    float pdf_wi = pdf_theta_i * pdf_phi_i / safe(cos(theta_i));
    Vec3f wi_local = to_cartesian(phi_i, theta_i, t_local);
    Vec3f wi_world = from_local(wi_local, normal);

    return {wi_world, pdf_wi};
  }

};

struct Microfacet {
    float roughness;
    // refraction indices for RGB channels
    Vec3f n1;
    Vec3f n2;

    /** Computes the Fresnel term for the microfacet material.
      \param[in] wi The light incident direction in local space.
      \return The Fresnel term.
    */
    Vec3f F(Vec3f wi) const {
        // =============================================================================================
        // TODO: Implement this function
        // check cosine
        Vec3f temp = (n1 - n2) / (n1 + n2);
        Vec3f R_0{temp[0] * temp[0], temp[1] * temp[1], temp[2] * temp[2]};

        float cos_theta = dot(wi, Vec3f{0.0f, 0.0f, 1.0f});
        return R_0 + (1 - R_0) * (float) (pow(1 - cos_theta, 5));
        // =============================================================================================
    }
    /** Computes the Beckmann normal distribution function for the microfacet material.
      \param[in] h The half vector in local space.
      \return The normal distribution function.
    */
    float D(Vec3f h) const {
        // =============================================================================================
        // TODO: Implement this function
        float cos_theta_h = dot(h, Vec3f{0.0f, 0.0f, 1.0f});

        float safe_cos_theta_h = cos_theta_h == 0 ? EPS : cos_theta_h;

        // if (std::isinf(tan(acos(cos_theta_h)))) {
        //   exit(0);
        // }
        // spdlog::info("{}", exp(-pow(tan(std::max(EPS, acos(cos_theta_h))), 2) / pow(roughness, 2)) / (M_PI * pow(roughness, 2) * pow(cos_theta_h, 4)));
        // spdlog::info("{}",exp(-pow(tan(std::max(EPS, acos(cos_theta_h))), 2) / pow(roughness, 2)));

        return exp(-pow(tan(acos(cos_theta_h)), 2) / pow(roughness, 2)) / (M_PI * pow(roughness, 2) * pow(safe_cos_theta_h, 4));
        // =============================================================================================
    }


    /** Computes the shadowing-masking function for the microfacet material.
      \param[in] wo The outgoing direction in local space.
      \param[in] wi The light incident direction in local space.
      \return The shadowing-masking value.
    */
    float G(Vec3f wo, Vec3f wi) const {
        // =============================================================================================
        // TODO: Implement this function
        float theta_o = acos(dot(wo, Vec3f{0.0f, 0.0f, 1.0f}));
        float theta_i = acos(dot(wi, Vec3f{0.0f, 0.0f, 1.0f}));

        // if (std::isinf(tan(theta_o)) || std::isinf(tan(theta_i))) {
        //   exit(0);
        // }

        float a_o = 1 / (roughness * tan(theta_o));
        float a_i = 1 / (roughness * tan(theta_i));

        float safe_a_o = a_o == 0 ? EPS : a_o;
        float safe_a_i = a_i == 0 ? EPS : a_i;


        float lambda_o = a_o < 1.6 ? (1 - 1.259 * a_o + 0.396 * pow(a_o, 2)) / (3.535 * safe_a_o + 2.181 * pow(safe_a_o, 2)) : 0;
        float lambda_i = a_i < 1.6 ? (1 - 1.259 * a_i + 0.396 * pow(a_i, 2)) / (3.535 * safe_a_i + 2.181 * pow(safe_a_i, 2)) : 0;

        float G_wo = 1 / (1 + lambda_o);
        float G_wi = 1 / (1 + lambda_i);

        return G_wo * G_wi;
        // =============================================================================================
    }
    /** Evaluates the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] wi_world The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The BRDF (fr) value.
    */
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        Vec3f wo = to_local(wo_world, normal);
        Vec3f wi = to_local(wi_world, normal);
        // =============================================================================================
        // TODO: Implement this function

        Vec3f wh = normalize(wo + wi);
        Vec3f n{0.0f, 0.0f, 1.0f};

        Vec3f f = F(wi);
        float d = D(wh);
        float g = G(wo, wi);
        
        // if ((4 * dot(wo, n) * dot(wi, n)) <= 0) {
        //   spdlog::info("{}", 4 * dot(wo, n) * dot(wi, n));
        //   exit(0);
        // }

        Vec3f fr = f * d * g / (4 * dot(wo, n) * dot(wi, n));
        // if (fr[0] < 0 || fr[1] < 0 || fr[2] < 0) {
        //   exit(0);
        // }

        return fr;
        // =============================================================================================
    }
    /** Computes the PDF for the microfacet material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        // =============================================================================================
        // TODO: Implement this function
        Vec3f wo = to_local(wo_world, normal);
        Vec3f wi = to_local(wi_world, normal);

        Vec3f wh = normalize(wo + wi);
        
        float phi_h = atan(wh.y / wh.x);

        float p_phi_h = 1 / (2 * M_PI);
        
        float theta_h = acos(wh.z / sqrt(length_squared(wh)));

        // if (theta_h == 0) {
        //   theta_h = EPS;
        // }

        float safeSin = theta_h == 0 ? EPS : sin(theta_h);
        float safeTan = cos(theta_h) == 0 ? tan(theta_h - EPS) : tan(theta_h);
        float safeCos = cos(theta_h) == 0 ? EPS : cos(theta_h);

        // float p_theta_h = 2 * safeSin / (pow(roughness, 2) * pow(posCos, 3)) * exp(-pow(safeTan, 2) / pow(roughness, 2));
        float p_theta_h = 2 * safeSin / (pow(roughness, 2) * pow(safeCos, 3)) * exp(-pow(safeTan, 2) / pow(roughness, 2));

        // if (std::isinf(p_theta_h) || p_theta_h <= 0) {
        //   spdlog::info("{}", theta_h);
        //   spdlog::info("{}", p_theta_h);
        //   // spdlog::info("{}", sqrt(length_squared(h)));
        //   // spdlog::info("{}", phi_h);
        //   // spdlog::info("{}", theta_h);
        //   exit(0);
        // }

        float p_omega_h = p_theta_h * p_phi_h / safeSin;

        float p_omega_i = p_omega_h / std::max(EPS, 4 * dot(wo, wh));
        
        return p_omega_i;
        // =============================================================================================
    }

    /** Samples the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
     std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u) const {
        // =============================================================================================
        // TODO: Implement this function
        Vec3f wo = to_local(wo_world, normal);

        float phi_h = 2 * M_PI * u[0];
        float theta_h = atan(sqrt(std::max((double) 0, -pow(roughness, 2) * log(1 - u[1])))); // check negative
        // float theta_h = atan(sqrt(-pow(roughness, 2) * log(1 - u[1]))); // check negative
        // spdlog::info("hi{}", phi_h);
        // spdlog::info("hi{}", theta_h);

        Vec3f wh = normalize(Vec3f{sin(theta_h) * cos(phi_h), sin(phi_h) * sin(theta_h), cos(theta_h)});

        // spdlog::info("w_h {}", w_h);
        // spdlog::info("wo_world {}", wo_world);

        Vec3f wi = normalize(mirror_reflect(-wo, wh));

        // spdlog::info("wi {}", wi);
        
        Vec3f wi_world = from_local(wi, normal);

        float pdf_wi = pdf(wo_world, wi_world, normal);

        return {wi_world, pdf_wi};
        // =============================================================================================
      }

};
}  // namespace muni
