#include "common.h"
#include "triangle.h"
#include "material.h"
#include "math_helpers.h"
#include <array>
#include <variant>

namespace muni { namespace BoxScene {

static const float light_x = 0.195f;
static const float light_y = -0.355f;
static const float light_z = 0.545f;
static const float light_len_x = 0.16f;
static const float light_len_y = 0.16f;
static const float inv_light_area = 1 / (light_len_x * light_len_y);
static const Vec3f light_color{50.0f, 50.0f, 50.0f};
static const Vec3f light_normal{0.0f, 0.0f, -1.0f};


const Cloth Cloth1{
    .A1 = Vec3f{0.06f, 0.24f, 0.3f},
    .A2 = Vec3f{0.06f, 0.24f, 0.3f},
    .gamma_s1 = 12 * M_PI / 180,
    .gamma_s2 = 12 * M_PI / 180,
    .gamma_v1 = 24 * M_PI / 180,
    .gamma_v2 = 24 * M_PI / 180,
    .k_d1 = 0.3,
    .k_d2 = 0.3,
    .eta3 = Vec3f{1.46f},
    .a1 = 0.33,
    .a2 = 0.33,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-25.0f, 25.0f},
    .t1_offset_lengths = std::vector<float>{1.0f},
    .t2_offsets = std::vector<float>{-25.0f, 25.0f},
    .t2_offset_lengths = std::vector<float>{1.0f},
    .N1 = 3,
    .N2 = 2
};


const Cloth Cloth2{
    .A1 = Vec3f{0.12f, 0.114f, 0.006f},
    .A2 = Vec3f{0.16f, 0.152f, 0.008f},
    .gamma_s1 = 5 * M_PI / 180,
    .gamma_s2 = 18 * M_PI / 180,
    .gamma_v1 = 10 * M_PI / 180,
    .gamma_v2 = 32 * M_PI / 180,
    .k_d1 = 0.2,
    .k_d2 = 0.3,
    .eta3 = Vec3f{1.345f},
    .a1 = 0.75,
    .a2 = 0.25,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-35.0f, -35.0f, 35.0f, 35.0f},
    .t1_offset_lengths = std::vector<float>{1.0f, 1.0f, 1.0f},
    .t2_offsets = std::vector<float>{0.0f, 0.0f},
    .t2_offset_lengths = std::vector<float>{1.0f},
    .N1 = 5,
    .N2 = 4
};

const Cloth Cloth3{
    .A1 = Vec3f{0.035, 0.01295f, 0.0105f},
    .A2 = Vec3f{0.2f, 0.074f, 0.06f},
    .gamma_s1 = 2.5 * M_PI / 180,
    .gamma_s2 = 30 * M_PI / 180,
    .gamma_v1 = 5 * M_PI / 180,
    .gamma_v2 = 60 * M_PI / 180,
    .k_d1 = 0.1,
    .k_d2 = 0.7,
    .eta3 = Vec3f{1.539f},
    .a1 = 0.9,
    .a2 = 0.1,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-32.0f, -32.0f, -18.0f, 0.0f, 0.0f, 18.0f, 32.0f, 32.0f},
    .t1_offset_lengths = std::vector<float>{1.33f, 0.66f, 2.0f, 2.0f, 2.0f, 0.66f, 1.33f},
    .t2_offsets = std::vector<float>{0.0f, 0.0f},
    .t2_offset_lengths = std::vector<float>{1.0f},
    .N1 = 6,
    .N2 = 2
};

const Cloth Cloth4{
    .A1 = Vec3f{0.035, 0.01295f, 0.0105f},
    .A2 = Vec3f{0.2f, 0.074f, 0.06f},
    .gamma_s1 = 2.5 * M_PI / 180,
    .gamma_s2 = 30 * M_PI / 180,
    .gamma_v1 = 5 * M_PI / 180,
    .gamma_v2 = 60 * M_PI / 180,
    .k_d1 = 0.1,
    .k_d2 = 0.7,
    .eta3 = Vec3f{1.539f},
    .a1 = 0.67,
    .a2 = 0.33,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-30.0f, -30.0f, 30.0f, 30.0f, -5.0f, -5.0f, 5.0f, 5.0f},
    .t1_offset_lengths = std::vector<float>{1.33f, 1.33f, 1.33f, 0.0f, 0.67f, 0.67f, 0.67f},
    .t2_offsets = std::vector<float>{0.0f, 0.0f},
    .t2_offset_lengths = std::vector<float>{1.0f},
    .N1 = 4,
    .N2 = 2
};

const Cloth Cloth5{
    .A1 = Vec3f{0.02f, 0.2f, 0.08f},
    .A2 = Vec3f{0.6f, 0.0f, 0.06f},
    .gamma_s1 = 4 * M_PI / 180,
    .gamma_s2 = 5 * M_PI / 180,
    .gamma_v1 = 8 * M_PI / 180,
    .gamma_v2 = 10 * M_PI / 180,
    .k_d1 = 0.1,
    .k_d2 = 0.1,
    .eta3 = Vec3f{1.345f},
    .a1 = 0.86,
    .a2 = 0.14,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-25.0f, -25.0f, 25.0f, 25.0f},
    .t1_offset_lengths = std::vector<float>{1.33f, 2.67f, 1.33f},
    .t2_offsets = std::vector<float>{0.0f, 0.0f},
    .t2_offset_lengths = std::vector<float>{1.0f},
    .N1 = 8,
    .N2 = 1
};

const Cloth Cloth6{
    .A1 = Vec3f{0.015f, 0.006f, 0.0f},
    .A2 = Vec3f{0.015f, 0.006f, 0.0f},
    .gamma_s1 = 6 * M_PI / 180,
    .gamma_s2 = 6 * M_PI / 180,
    .gamma_v1 = 12 * M_PI / 180,
    .gamma_v2 = 12 * M_PI / 180,
    .k_d1 = 0.1,
    .k_d2 = 0.1,
    .eta3 = Vec3f{1.46f},
    .a1 = 0.5,
    .a2 = 0.5,
    .t1 = Vec3f{1.0f, 0.0f, 0.0f},
    .t2 = Vec3f{0.0f, 1.0f, 0.0f},
    .t1_offsets = std::vector<float>{-90.0f, 50.0f},
    .t1_offset_lengths = std::vector<float>{1.0f},
    .t2_offsets = std::vector<float>{-90.0f, -55.0f, 55.0f, 90.0f},
    .t2_offset_lengths = std::vector<float>{0.5f, EPS, 0.5f},
    .N1 = 5,
    .N2 = 4
};





const Microfacet Gold{.roughness = 0.0005f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{0.2177f, 0.42659f,1.2425f } };
const Microfacet Iron{.roughness = 0.02f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{2.8851f, 2.95f,2.6f } };
static const std::array<std::variant<Lambertian, Microfacet, Cloth>, 13> materials = {
    // Back
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Bottom
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Left
    Lambertian{.albedo = Vec3f{0.0f, 0.2117f, 0.3765f}},
    // Right
    Lambertian{.albedo = Vec3f{0.996f, 0.7373f, 0.0667f}},
    // Top
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Bunny
    Iron,
    Gold,
    Cloth1,
    Cloth2,
    Cloth3,
    Cloth4,
    Cloth5,
    Cloth6
};

static std::vector<Triangle> triangles = {
    // Light
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y, light_z},
             .v2 = Vec3f{light_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y + light_len_y, light_z},
             .v2 = Vec3f{light_x + light_len_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    // Back
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    // Bottom
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    // Left
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    // Right
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    // Top
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    };
}}  // namespace muni::BoxScene
