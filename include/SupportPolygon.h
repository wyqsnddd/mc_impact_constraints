#pragma once

#include <mc_rbdyn/Robot.h>
#include <mc_rbdyn/SCHAddon.h>
#include <mc_rbdyn/Surface.h>
#include <mc_rbdyn/surface_hull.h>

#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullPoints.h>
#include <libqhullcpp/QhullVertexSet.h>

namespace mc_impact
{

class SupportPolygon
{
public:
  SupportPolygon(const mc_rbdyn::Robot & robot) : robot_(robot)
  {
    std::cout << "Support polygon calculator is constructed." << std::endl;
  }
  ~SupportPolygon()
  {
    free(G_zmp_);
    free(h_zmp_);
    free(supportPolygonHullPtr_);
  }

  void addContactSurface(const std::string &);
  void removeContactSurface(const std::string &);
  std::vector<std::string> existingSurfaces() const;
  /*
  inline const std::vector< sva::PTransformd > & verticies() const
  {
    return verticies_;
  }
  */
  void update();
  bool readMatricies(Eigen::MatrixXd & G, Eigen::VectorXd & h);
  inline const sch::S_Polyhedron * getSupportPolygon()
  {
    return supportPolygonHullPtr_;
  }

private:
  const mc_rbdyn::Robot & robot_;
  std::map<std::string, const mc_rbdyn::Surface &> surfaceContainer_;
  // std::vector< sva::PTransformd > verticies_;
  Eigen::MatrixXd * G_zmp_;
  Eigen::VectorXd * h_zmp_;
  sch::S_Polyhedron * supportPolygonHullPtr_;
  sch::S_Polyhedron * createPoly_(const orgQhull::Qhull & qhull);
};

} // namespace mc_impact
