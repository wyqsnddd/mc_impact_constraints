#pragma once

#include <mc_rtc/GUIState.h>
#include <mc_solver/InequalityConstraint.h>
#include <mc_solver/QPSolver.h>

#include "ContactWrenchMatrixToLambdaMatrix.h"


namespace mc_impact
{

/** This constraint enforce a minimum wrench for the QP forces.
 *
 * This is meant as a simple for ContactWrenchMatrixToLambdaMatrix.
 *
 * It enforces \f$-F_{QP} < -F_{user}\f$ i.e. \f$ F_{QP} > F_{user} \f$ where \f$F_{user}\f$ can be changed by the user
 */
struct MinimumContactWrench : public mc_solver::InequalityConstraint
{
  MinimumContactWrench(const mc_solver::QPSolver & solver,
                       const mc_rbdyn::Contact & contact,
                       const sva::ForceVecd & fUser,
                       std::shared_ptr<mc_rtc::gui::StateBuilder> gui)
  : InequalityConstraint(contact.contactId(solver.robots()))
  {
    auto cid = contact.contactId(solver.robots());
    ContactWrenchMatrixToLambdaMatrix transformer(solver, cid);
    sName_ = contact.r1Surface()->name();
    A_ = transformer.transform(-Eigen::MatrixXd::Identity(6, 6));
    fUser_ = -fUser.vector();
    gui_ = gui;
    if(gui_)
    {
      gui_->addElement({"Constraints", "MinimumContactWrench", sName_},
                       mc_rtc::gui::ArrayInput("Minimum wrench", {"cx", "cy", "cz", "fx", "fy", "fz"}, [this]() -> Eigen::VectorXd { return -fUser_; }, [this](const Eigen::VectorXd & v) { fUser_ = -v; }));
    }
  }

  ~MinimumContactWrench() override
  {
    if(gui_)
    {
      gui_->removeCategory({"Constraints", "MinimumContactWrench", sName_});
    }
  }

  int maxInEq() const override { return 6; }

  std::string nameInEq() const override { return "MinimumContactWrench"; }

  const Eigen::VectorXd & bInEq() const override { return fUser_; }

  const Eigen::MatrixXd & A() const override { return A_; }

  void computeAb() override
  {
  }
private:
  Eigen::VectorXd fUser_;
  Eigen::MatrixXd A_;
  std::shared_ptr<mc_rtc::gui::StateBuilder> gui_;
  std::string sName_;
};

}
