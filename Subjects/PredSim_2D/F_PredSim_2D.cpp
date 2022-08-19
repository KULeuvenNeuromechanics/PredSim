/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
#include "SimTKcommon/internal/recorder.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace OpenSim;

/*  The function F describes the OpenSim model and, implicitly, the skeleton
    dynamics. F takes as inputs joint positions and velocities (states x),
    joint accelerations (controls u), and returns the joint torques as well as
    several variables for use in the optimal control problems. F is templatized
    using type T. F(x,u)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 10;        // # degrees of freedom
constexpr int NX = 2*ndof;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof + 4;    // # residual torques + # GRF

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// OpenSim and Simbody use different indices for the states/controls when the
// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
// and legs with pelvis as origin).
// The two following functions allow getting the indices from one reference
// system to the other. These functions are inspired from
// createSystemYIndexMap() in Moco.
// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
// reference system. Note that we only care about the order here so we divide
// by 2 because the states include both Qs and Qdots.
SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
    auto s = model.getWorkingState();
    const auto svNames = model.getStateVariableNames();
    SimTK::Array_<int> idxOSInSimbody(s.getNQ());
    s.updQ() = 0;
    for (int iy = 0; iy < s.getNQ(); ++iy) {
        s.updQ()[iy] = SimTK::NaN;
        const auto svValues = model.getStateVariableValues(s);
        for (int isv = 0; isv < svNames.size(); ++isv) {
            if (SimTK::isNaN(svValues[isv])) {
                s.updQ()[iy] = 0;
                idxOSInSimbody[iy] = isv/2;
                break;
            }
        }
    }
    return idxOSInSimbody;
}
// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
// reference system.
SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
    auto idxOSInSimbody = getIndicesOSInSimbody(model);
    auto s = model.getWorkingState();
    SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
			if (idxOSInSimbody[iyy] == iy) {
				idxSimbodyInOS[iy] = iyy;
				break;
			}
		}
	}
    return idxSimbodyInOS;
}

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    // OpenSim model: create components
    /// Model
    OpenSim::Model* model;
    /// Bodies
    OpenSim::Body* pelvis;
    OpenSim::Body* femur_r;
    OpenSim::Body* femur_l;
    OpenSim::Body* tibia_r;
    OpenSim::Body* tibia_l;
    OpenSim::Body* talus_r;
    OpenSim::Body* talus_l;
    OpenSim::Body* calcn_r;
    OpenSim::Body* calcn_l;
    OpenSim::Body* toes_r;
    OpenSim::Body* toes_l;
    OpenSim::Body* torso;
    /// Joints
    OpenSim::PlanarJoint* ground_pelvis;
    OpenSim::CustomJoint* hip_r;
    OpenSim::CustomJoint* hip_l;
    OpenSim::CustomJoint* knee_r;
    OpenSim::CustomJoint* knee_l;
    OpenSim::PinJoint* ankle_r;
    OpenSim::PinJoint* ankle_l;
    OpenSim::WeldJoint* subtalar_r;
    OpenSim::WeldJoint* subtalar_l;
    OpenSim::WeldJoint* mtp_r;
    OpenSim::WeldJoint* mtp_l;
    OpenSim::PinJoint* back;
    /// Contact elements
    OpenSim::HuntCrossleyForce_smooth* HC_heel_r;
    OpenSim::HuntCrossleyForce_smooth* HC_front_r;
    OpenSim::HuntCrossleyForce_smooth* HC_heel_l;
    OpenSim::HuntCrossleyForce_smooth* HC_front_l;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 9.7143336091724, Vec3(-0.0682778, 0, 0), Inertia(0.0814928846050306, 0.0814928846050306, 0.0445427591530667, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 7.67231915023828, Vec3(0, -0.170467, 0), Inertia(0.111055472890139, 0.0291116288158616, 0.117110028170931, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 7.67231915023828, Vec3(0, -0.170467, 0), Inertia(0.111055472890139, 0.0291116288158616, 0.117110028170931, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 3.05815503574821, Vec3(0, -0.180489, 0), Inertia(0.0388526996597354, 0.00393152317985418, 0.0393923204883429, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 3.05815503574821, Vec3(0, -0.180489, 0), Inertia(0.0388526996597354, 0.00393152317985418, 0.0393923204883429, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.082485638186061, Vec3(0), Inertia(0.000688967700910182, 0.000688967700910182, 0.000688967700910182, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.082485638186061, Vec3(0), Inertia(0.000688967700910182, 0.000688967700910182, 0.000688967700910182, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 1.03107047732576, Vec3(0.0913924, 0.0274177, 0), Inertia(0.000964554781274254, 0.00268697403354971, 0.00282476757373175, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 1.03107047732576, Vec3(0.0913924, 0.0274177, 0), Inertia(0.000964554781274254, 0.00268697403354971, 0.00282476757373175, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.178663892311008, Vec3(0.0316218, 0.00548355, 0.0159937), Inertia(6.88967700910182e-005, 0.000137793540182036, 6.88967700910182e-005, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.178663892311008, Vec3(0.0316218, 0.00548355, -0.0159937), Inertia(6.88967700910182e-005, 0.000137793540182036, 6.88967700910182e-005, 0, 0, 0));
    torso = new OpenSim::Body("torso", 28.240278003209, Vec3(-0.0289722, 0.309037, 0), Inertia(1.14043571182129, 0.593400919285897, 1.14043571182129, 0, 0, 0));
    /// Joints
    /// Hip_l transform
    SpatialTransform st_hip_l;
    st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
    st_hip_l[2].setFunction(new LinearFunction());
    st_hip_l[2].setAxis(Vec3(0, 0, 1));
    /// Hip_r transform
    SpatialTransform st_hip_r;
    st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
    st_hip_r[2].setFunction(new LinearFunction());
    st_hip_r[2].setAxis(Vec3(0, 0, 1));
    /// Knee_l transform
    SpatialTransform st_knee_l;
    st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    st_knee_l[2].setFunction(new LinearFunction());
    st_knee_l[2].setAxis(Vec3(0, 0, 1));
    /// Knee_r transform
    SpatialTransform st_knee_r;
    st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[2].setFunction(new LinearFunction());
    st_knee_r[2].setAxis(Vec3(0, 0, 1));
    /// Joint specifications
    ground_pelvis = new PlanarJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0));
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.0682778001711179, -0.0638353973311301, -0.0823306940058688), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.0682778001711179, -0.0638353973311301, 0.0823306940058688), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(-0.00451221232146798, -0.396907245921447, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.00451221232146798, -0.396907245921447, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    ankle_l = new PinJoint("ankle_l", *tibia_l, Vec3(0, -0.415694825374905, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0));
    ankle_r = new PinJoint("ankle_r", *tibia_r, Vec3(0, -0.415694825374905, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0));
    subtalar_l = new WeldJoint("subtalar_l", *talus_l, Vec3(-0.0445720919117321, -0.0383391276542374, -0.00723828107321956), Vec3(0), *calcn_l, Vec3(0), Vec3(0));
    subtalar_r = new WeldJoint("subtalar_r", *talus_r, Vec3(-0.0445720919117321, -0.0383391276542374, 0.00723828107321956), Vec3(0), *calcn_r, Vec3(0), Vec3(0));
    mtp_l = new WeldJoint("mtp_l", *calcn_l, Vec3(0.163409678774199, -0.00182784875586352, -0.000987038328166303), Vec3(0), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new WeldJoint("mtp_r", *calcn_r, Vec3(0.163409678774199, -0.00182784875586352, 0.000987038328166303), Vec3(0), *toes_r, Vec3(0), Vec3(0));
    back = new PinJoint("back", *pelvis, Vec3(-0.0972499926058214, 0.0787077894476112, 0), Vec3(0), *torso, Vec3(0), Vec3(0));
    /// Add bodies and joints to model
    model->addBody(pelvis);		    model->addJoint(ground_pelvis);
    model->addBody(femur_l);		model->addJoint(hip_l);
    model->addBody(femur_r);		model->addJoint(hip_r);
    model->addBody(tibia_l);		model->addJoint(knee_l);
    model->addBody(tibia_r);		model->addJoint(knee_r);
    model->addBody(talus_l);		model->addJoint(ankle_l);
    model->addBody(talus_r);		model->addJoint(ankle_r);
    model->addBody(calcn_l);		model->addJoint(subtalar_l);
    model->addBody(calcn_r);		model->addJoint(subtalar_r);
    model->addBody(toes_l);		    model->addJoint(mtp_l);
    model->addBody(toes_r);		    model->addJoint(mtp_r);
    model->addBody(torso);          model->addJoint(back);
    /// Contact elements
    /// Parameters
    osim_double_adouble radiusSphere_heel = 0.035;
    osim_double_adouble radiussphere_front = 0.015;
    osim_double_adouble stiffness_heel = 3067776;
    osim_double_adouble stiffness_front = 3067776;
    osim_double_adouble dissipation = 2.0;
    osim_double_adouble staticFriction = 0.8;
    osim_double_adouble dynamicFriction = 0.8;
    osim_double_adouble viscousFriction = 0.5;
    osim_double_adouble transitionVelocity = 0.2;
    Vec3 locSphere_heel_l = Vec3(0.031307527581931796, 0.010435842527310599, 0);
    Vec3 locsphere_front_l = Vec3(0.1774093229642802, -0.015653763790965898, -0.005217921263655299);
    Vec3 locSphere_heel_r = Vec3(0.031307527581931796, 0.010435842527310599, 0);
    Vec3 locsphere_front_r = Vec3(0.1774093229642802, -0.015653763790965898, 0.005217921263655299);
    Vec3 normal = Vec3(0, 1, 0);
    osim_double_adouble offset = 0;
    /// Right foot contact shere specifications
    HC_heel_r = new HuntCrossleyForce_smooth("sphere_heel_r", "calcn_r", locSphere_heel_r, radiusSphere_heel,
        stiffness_heel, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_front_r = new HuntCrossleyForce_smooth("sphere_front_r", "calcn_r", locsphere_front_r, radiussphere_front,
        stiffness_front, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add right foot contact spheres to model
    model->addComponent(HC_heel_r);
    HC_heel_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_front_r);
    HC_front_r->connectSocket_body_sphere(*calcn_r);
    /// Left foot contact shere specifications
    HC_heel_l = new HuntCrossleyForce_smooth("sphere_heel_l", "calcn_l", locSphere_heel_l, radiusSphere_heel,
        stiffness_heel, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_front_l = new HuntCrossleyForce_smooth("sphere_front_l", "calcn_l", locsphere_front_l, radiussphere_front,
        stiffness_front, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add left foot contact spheres to model
    model->addComponent(HC_heel_l);
    HC_heel_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_front_l);
    HC_front_l->connectSocket_body_sphere(*calcn_l);

    // Initialize system and state.
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);

    // States and controls
    T ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    QsUs.setToZero();
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Compute residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndof);
    appliedMobilityForces.setToZero();
    /// appliedBodyForces (# bodies + ground)
    Vector_<SpatialVec> appliedBodyForces;
    int nbodies = model->getBodySet().getSize() + 1;
    appliedBodyForces.resize(nbodies);
    appliedBodyForces.setToZero();
    /// Set gravity
    Vec3 gravity(0);
    gravity[1] = -9.81;
    /// Add weights to appliedBodyForces
    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
        model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
    }
    /// Add contact forces to appliedBodyForces
    /// Right foot
    Array<osim_double_adouble> Force_values_heel_r = HC_heel_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_front_r = HC_front_r->getRecordValues(*state);
    SpatialVec GRF_heel_r;
    GRF_heel_r[0] = Vec3(Force_values_heel_r[9], Force_values_heel_r[10], Force_values_heel_r[11]);
    GRF_heel_r[1] = Vec3(Force_values_heel_r[6], Force_values_heel_r[7], Force_values_heel_r[8]);
    SpatialVec GRF_front_r;
    GRF_front_r[0] = Vec3(Force_values_front_r[9], Force_values_front_r[10], Force_values_front_r[11]);
    GRF_front_r[1] = Vec3(Force_values_front_r[6], Force_values_front_r[7], Force_values_front_r[8]);
    int nfoot_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
    appliedBodyForces[nfoot_r] = appliedBodyForces[nfoot_r] + GRF_heel_r + GRF_front_r;
    /// Left foot
    Array<osim_double_adouble> Force_values_heel_l = HC_heel_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_front_l = HC_front_l->getRecordValues(*state);
    SpatialVec GRF_heel_l;
    GRF_heel_l[0] = Vec3(Force_values_heel_l[9], Force_values_heel_l[10], Force_values_heel_l[11]);
    GRF_heel_l[1] = Vec3(Force_values_heel_l[6], Force_values_heel_l[7], Force_values_heel_l[8]);
    SpatialVec GRF_front_l;
    GRF_front_l[0] = Vec3(Force_values_front_l[9], Force_values_front_l[10], Force_values_front_l[11]);
    GRF_front_l[1] = Vec3(Force_values_front_l[6], Force_values_front_l[7], Force_values_front_l[8]);
    int nfoot_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    appliedBodyForces[nfoot_l] = appliedBodyForces[nfoot_l] + GRF_heel_l + GRF_front_l;
    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
        *state, appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Ground reaction forces
    SpatialVec GRF_r = GRF_heel_r + GRF_front_r;
    SpatialVec GRF_l = GRF_heel_l + GRF_front_l;

    // Extract results
    /// Residual forces
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < NU; ++i) res[0][i] =
            value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    for (int i = 0; i < 2; ++i) {
            res[0][i + ndof] = value<T>(GRF_r[1][i]); // GRF_r (x and y)
        }
        for (int i = 0; i < 2; ++i) {
            res[0][i + ndof + 2] = value<T>(GRF_l[1][i]); // GRF_l (x and y)
        }

    return 0;

}

/* In main(), the Recorder is used to save the expression graph of function F.
This expression graph is saved as a MATLAB function named foo.m. From this
function, a c-code can be generated via CasADi and then compiled as a dll. This
dll is then imported in MATLAB as an external function. With this workflow,
CasADi can use algorithmic differentiation to differentiate the function F.
*/
int main() {

    Recorder x[NX];
    Recorder u[NU];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;
}

