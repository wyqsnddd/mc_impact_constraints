#include "SupportPolygon.h"

namespace mc_impact
{

void SupportPolygon::addContactSurface(const std::string & surfaceName)
{
  // Find the surface
  if(!robot_.hasSurface(surfaceName))
    throw std::runtime_error(std::string("supportPolygon failed to find the contact surface: ") + surfaceName);
  surfaceContainer_.insert(std::pair<std::string, const mc_rbdyn::Surface &>(surfaceName, robot_.surface(surfaceName)));
}

void SupportPolygon::removeContactSurface(const std::string & surfaceName)
{
  const auto ee = surfaceContainer_.find(surfaceName);

  if(ee == surfaceContainer_.end())
  {
    throw std::runtime_error(std::string("Surface name: ") + surfaceName + std::string(" is not found. "));
  }
  else
  {
    surfaceContainer_.erase(ee);
  }
}

std::vector<std::string> SupportPolygon::existingSurfaces() const
{
  std::vector<std::string> * nameVectors;
  for(auto it = surfaceContainer_.begin(); it != surfaceContainer_.end(); ++it)
  {
    nameVectors->push_back(it->first);
  }
  return *nameVectors;
}

void SupportPolygon::update()
{

  // sch::S_Polyhedron * poly = new sch::S_Polyhedron();
  // auto & poly_algo = *(poly->getPolyhedronAlgorithm());

  int points_size = 0;
  // Go through the surfaces
  for(auto surfaceIdx = surfaceContainer_.begin(); surfaceIdx != surfaceContainer_.end(); ++surfaceIdx)
  {
    points_size += surfaceIdx->second.points().size();
  }
  // Go through the points: update the verticies_
  std::vector<double> points_in;
  points_in.reserve(points_size * 3);

  for(auto surfaceIdx = surfaceContainer_.begin(); surfaceIdx != surfaceContainer_.end(); ++surfaceIdx)
  {
    for(const auto & p : surfaceIdx->second.points())
    {
      const auto & t = p.translation();
      points_in.push_back(t.x());
      points_in.push_back(t.y());
      points_in.push_back(t.z());
    }
  } // end of for points
  orgQhull::Qhull qhull;
  qhull.runQhull("", 3, points_size, points_in.data(), "Qt");
  // auto points = qhull.points();
  supportPolygonHullPtr_ = createPoly_(qhull);

  int vertexNumer = qhull.vertexCount();

  G_zmp_ = new Eigen::MatrixXd(vertexNumer, 2);
  h_zmp_ = new Eigen::VectorXd(vertexNumer);
  // Eigen::MatrixXd G(vertexNumer, 2);
  // Eigen::VectorXd h(vertexNumer);
  G_zmp_->setOnes();
  h_zmp_->setOnes();

  auto tempVertexList = qhull.vertexList();
  int vNumber = 0;
  Eigen::Vector3d center;
  center << qhull.feasiblePoint()[0], qhull.feasiblePoint()[1], qhull.feasiblePoint()[2];

  // poly_algo.vertexes_.reserve(points.size());
  for(auto ii = tempVertexList.begin(); ii != tempVertexList.end(); ++ii, ++vNumber)
  {
    Eigen::Vector3d point_one, point_two;
    point_one << ii->point().coordinates()[0], ii->point().coordinates()[1], ii->point().coordinates()[2];

    if((ii + 1) == tempVertexList.end())
    {
      auto jj = tempVertexList.begin();
      point_two << jj->point().coordinates()[0], jj->point().coordinates()[1], jj->point().coordinates()[2];
    }
    else
    {
      point_two << (ii + 1)->point().coordinates()[0], (ii + 1)->point().coordinates()[1],
          (ii + 1)->point().coordinates()[2];
    }

    Eigen::Vector3d difference = point_two - point_one;
    difference.normalize();
    double slope = difference.y() / difference.x();

    (*G_zmp_)(vNumber, 0) = -slope;
    (*h_zmp_)(vNumber) = -slope * point_one.x() + point_one.y();

    int lineSign = 1;
    if(!(center.y() - slope * center.x() <= -slope * ((*h_zmp_)(vNumber))))
    {
      lineSign = -1;
    }
    // Correct the sign with the centeroid point
    (*G_zmp_).block(vNumber, 0, 1, 2) *= lineSign;
    (*h_zmp_)(vNumber) *= lineSign;
  } // iterate for each vertex
}

bool SupportPolygon::readMatricies(Eigen::MatrixXd & G, Eigen::VectorXd & h)
{
  if((G_zmp_ == nullptr) || (h_zmp_ == nullptr))
    return false;
  else
  {
    G = *G_zmp_;
    h = *h_zmp_;
    return true;
  }
}

sch::S_Polyhedron * SupportPolygon::createPoly_(const orgQhull::Qhull & qhull)
{
  sch::S_Polyhedron * poly = new sch::S_Polyhedron();
  auto & poly_algo = *(poly->getPolyhedronAlgorithm());

  auto points = qhull.points();
  poly_algo.vertexes_.reserve(points.size());
  for(const auto & p : points)
  {
    auto v = new sch::S_PolyhedronVertex();
    v->setCoordinates(p.coordinates()[0], p.coordinates()[1], p.coordinates()[2]);
    v->setNumber(p.id());
    poly_algo.vertexes_.push_back(v);
  }
  auto facets = qhull.facetList();
  poly_algo.triangles_.reserve(facets.size());
  for(const auto & f : facets)
  {
    if(!f.isGood())
    {
      continue;
    }
    sch::PolyhedronTriangle t;
    t.normal.Set(f.hyperplane().coordinates());
    t.normal.normalize();
    t.a = f.vertices()[0].point().id();
    t.b = f.vertices()[1].point().id();
    t.c = f.vertices()[2].point().id();
    auto addNeighbors = [](std::vector<sch::S_PolyhedronVertex *> & vertexes_, unsigned int a, unsigned int b,
                           unsigned int c) {
      vertexes_[a]->addNeighbor(vertexes_[b]);
      vertexes_[a]->addNeighbor(vertexes_[c]);
    };
    addNeighbors(poly_algo.vertexes_, t.a, t.b, t.c);
    addNeighbors(poly_algo.vertexes_, t.b, t.a, t.c);
    addNeighbors(poly_algo.vertexes_, t.c, t.a, t.b);
    poly_algo.triangles_.push_back(t);
  }

  for(const auto & v : poly_algo.vertexes_)
  {
    v->updateFastArrays();
  }
  poly_algo.deleteVertexesWithoutNeighbors();

  return poly;
}
} // namespace mc_impact
