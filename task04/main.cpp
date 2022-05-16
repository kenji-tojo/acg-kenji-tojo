/*
 * Copyright (c) 2019 Nobuyuki Umetani
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include <vector>
#include <algorithm>
#if defined(_WIN32) // windows
#  define NOMINMAX   // to remove min,max macro
#  include <windows.h>  // this should come before glfw3.h
#endif
#define GL_SILENCE_DEPRECATION
#include <GLFW/glfw3.h>

#include "delfem2/srch_trimesh3_class.h"
#include "delfem2/srch_bruteforce.h"
#include "delfem2/srch_bv3_sphere.h"
#include "delfem2/srch_bvh.h"
#include "delfem2/msh_affine_transformation.h"
#include "delfem2/msh_io_ply.h"
#include "delfem2/msh_normal.h"
#include "delfem2/mat4.h"
#include "delfem2/sampling.h"
#include "delfem2/opengl/tex.h"
#include "delfem2/glfw/viewer3.h"
#include "delfem2/glfw/util.h"

namespace dfm2 = delfem2;

// ----------------------------------------

/*! function to compute "weight" and "direction"
 * @param[out] dir
 * @param Xi random number generator
 * @param nrm normal of surface
 * @return weight
 */
double SamplingHemisphere(
    double dir[3],
    std::array<unsigned short, 3> &Xi,
    const double nrm[3])  //
{
  // below implement code to sample hemisphere with cosine weighted probabilistic distribution
  // hint1: use polar coordinate (longitude and latitude).
  // hint2: generate two float values using "dfm2::MyERand48<double>(Xi)". One will be longitude and another will be latitude
  // hint3: for longitude use inverse sampling method to achieve cosine weighted sample.
  // hint4: first assume z is the up in the polar coordinate, then rotate the sampled direction such that "z" will be up.
  // write some codes below (5-10 lines)
  
  // define local frame around nrm
  double local_x[3]{ 0 };
  if (fabs(nrm[0]) < 1e-5) {
    local_x[0] = 1.0;
  } else {
    const auto len = sqrt(nrm[0]*nrm[0]+nrm[1]*nrm[1]);
    local_x[0] =  nrm[1]/len;
    local_x[1] = -nrm[0]/len;
  }
  double local_y[3]{ // local_y = cross(nrm,local_x)
    /*local_x[2]==0*/ - nrm[2]*local_x[1],
    nrm[2]*local_x[0],
    nrm[0]*local_x[1] - nrm[1]*local_x[0]
  };
  const auto cos_theta = dfm2::MyERand48<double>(Xi);
  const auto sin_theta = sqrt(fmax(0,1-cos_theta*cos_theta));
  const auto phi = dfm2::MyERand48<double>(Xi) * M_PI * 2.0;
  const auto cos_phi = cos(phi);
  const auto sin_phi = sin(phi);
  for (int i=0;i<3;i++) {
    dir[i] = cos_theta*nrm[i] + sin_theta*cos_phi*local_x[i] + sin_theta*sin_phi*local_y[i];
  }

  return 1;
}

double SampleAmbientOcclusion(
    std::array<unsigned short, 3> &Xi,
    const dfm2::CVec3d &src1,
    const dfm2::CVec3d &dir1,
    const std::vector<double> &vec_xyz,
    const std::vector<unsigned int> &vec_tri,
    const std::vector<dfm2::CNodeBVH2> &bvh_nodes,
    const std::vector<dfm2::CBV3_Sphere<double>> &bvh_volumes) {
  dfm2::PointOnSurfaceMesh<double> pos_mesh;
  bool is_hit = Intersection_Ray3_Tri3_Bvh(
      pos_mesh,
      src1, dir1, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
  if (!is_hit) { return 0; } // the ray from pixel doesn't hit the mesh
  // ---------------
  unsigned int itri = pos_mesh.itri;  // the triangle hit by the ray
  assert(itri < vec_tri.size() / 3);
  dfm2::CVec3d nrm_tri = dfm2::Normal_TriInMeshTri3(itri, vec_xyz.data(), vec_tri.data());
  nrm_tri.normalize();
  dfm2::CVec3d src2 = pos_mesh.PositionOnMeshTri3(vec_xyz, vec_tri);
  src2 += nrm_tri * 1.0e-3;
  dfm2::CVec3d dir2;
  const double weight0 = ::SamplingHemisphere(
      dir2.p,
      Xi,
      nrm_tri.normalized().data());
  // check if the ray from the triangle hit the mesh
  dfm2::PointOnSurfaceMesh<double> pos_mesh2;
  bool is_hit2 = Intersection_Ray3_Tri3_Bvh(  // check if ray (src2,dir2) hit the triangle mesh
      pos_mesh2,
      src2, dir2, vec_xyz, vec_tri, bvh_nodes, bvh_volumes);
  if (!is_hit2) { return weight0; }
  return 0;
}

int main() {
  std::vector<double> vtx_xyz; // 3d points
  std::vector<unsigned int> tri_vtx;
  { // load input mesh
    delfem2::Read_Ply(
        vtx_xyz, tri_vtx,
        std::string(PATH_SRC_DIR) + "/../assets/bunny_2k.ply");
    dfm2::Rotate_Points3(vtx_xyz, -M_PI*0.5, 0.0, 0.0);
    dfm2::Normalize_Points3(vtx_xyz, 2.5);
  }
  // spatial hash data structure
  std::vector<dfm2::CNodeBVH2> bvh_nodes;
  std::vector<dfm2::CBV3_Sphere<double>> bvh_volumes;
  delfem2::ConstructBVHTriangleMeshMortonCode(
      bvh_nodes, bvh_volumes,
      vtx_xyz, tri_vtx);
  const unsigned int nw = 256;
  const unsigned int nh = 256;
  // above: constant
  // -------------------------------
  // below: changing during execution
  std::vector<float> afRGB(nw * nh * 3, 0.f);
  unsigned int isample = 0;
  dfm2::CMat4d mMVPd_inv;
  auto render = [&](int iw, int ih) {
    std::array<unsigned short, 3> Xi = {
        (unsigned short) (ih * ih),
        (unsigned short) (iw * iw),
        (unsigned short) (isample * isample)};
    const std::pair<dfm2::CVec3d,dfm2::CVec3d> ray = dfm2::RayFromInverseMvpMatrix(
        mMVPd_inv.data(), iw, ih, nw, nh );
    double rd0 = SampleAmbientOcclusion(
        Xi,
        ray.first, ray.second, vtx_xyz, tri_vtx, bvh_nodes, bvh_volumes);
    dfm2::CVec3d r_ave(rd0, rd0, rd0);
    {
      float *ptr = afRGB.data() + (ih * nw + iw) * 3;
      const auto isamplef = static_cast<float>(isample);
      ptr[0] = (isamplef * ptr[0] + r_ave[0]) / (isamplef + 1.f);
      ptr[1] = (isamplef * ptr[1] + r_ave[1]) / (isamplef + 1.f);
      ptr[2] = (isamplef * ptr[2] + r_ave[2]) / (isamplef + 1.f);
    }
  };

  dfm2::opengl::CTexRGB_Rect2D tex;
  {
    tex.width = nw;
    tex.height = nh;
    tex.channels = 3;
    tex.pixel_color.resize(tex.width * tex.height * tex.channels);
  }

  dfm2::glfw::CViewer3 viewer(2.f);
  viewer.window_title = "task04";
  viewer.width = 400;
  viewer.height = 400;
  viewer.camerachange_callbacks.emplace_back( // reset when camera moves
      [&afRGB, &isample] {
        std::fill(afRGB.begin(), afRGB.end(), 0.0);
        isample = 0;
      }
  );
  // --------------
  // start OpenGL

  if (!glfwInit()) { exit(EXIT_FAILURE); }
  ::glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
  ::glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  viewer.OpenWindow();
  tex.InitGL();

  while (!glfwWindowShouldClose(viewer.window)) {
    { // inverse of Homography matrix
      const dfm2::CMat4f mP = viewer.GetProjectionMatrix();
      const dfm2::CMat4f mMV = viewer.GetModelViewMatrix();
      const dfm2::CMat4d mMVP = (mP * mMV).cast<double>();
      mMVPd_inv = dfm2::Inverse_Mat4(mMVP.data());
    }
    /*
    for(unsigned int iw=0;iw<nw;++iw){
      for(unsigned int ih=0;ih<nh;++ih){
        render(iw,ih);
      }
    }
     */
    dfm2::parallel_for(nw, nh, render);
    isample++;
    for (unsigned int ih = 0; ih < tex.height; ++ih) {
      for (unsigned int iw = 0; iw < tex.width; ++iw) {
        for (int ic = 0; ic < 3; ++ic) {
          float fc = afRGB[(ih * tex.width + iw) * 3 + ic];
          fc = (fc > 1.f) ? 1.f : fc;
          fc = (fc < 0.f) ? 0.f : fc;
          int ifc = int(fc * 255.f + .5f);
          tex.pixel_color[(ih * tex.width + iw) * 3 + ic] = ifc;
        }
      }
    }
    tex.InitGL();
    ::glfwMakeContextCurrent(viewer.window);
    ::glClearColor(0.8, 1.0, 1.0, 1.0);
    ::glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D , tex.id_tex);
    tex.Draw_oldGL();
    viewer.SwapBuffers();
    glfwPollEvents();
  }

  glfwDestroyWindow(viewer.window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}


