/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/PythonForce.h"
#include "openmm/Platform.h"
#include "openmm/VerletIntegrator.h"
#include "sfmt/SFMT.h"
#include <cstdlib>
#include <iostream>
#include <thread>

using namespace OpenMM;
using namespace std;

class ScopedEnvironmentVariable {
public:
    ScopedEnvironmentVariable(const string& name, const string& value) : name(name), hadValue(false) {
        const char* currentValue = getenv(name.c_str());
        if (currentValue != nullptr) {
            hadValue = true;
            previousValue = currentValue;
        }
#ifdef WIN32
        _putenv_s(name.c_str(), value.c_str());
#else
        setenv(name.c_str(), value.c_str(), 1);
#endif
    }
    ~ScopedEnvironmentVariable() {
        if (hadValue) {
#ifdef WIN32
            _putenv_s(name.c_str(), previousValue.c_str());
#else
            setenv(name.c_str(), previousValue.c_str(), 1);
#endif
        }
        else {
#ifdef WIN32
            _putenv_s(name.c_str(), "");
#else
            unsetenv(name.c_str());
#endif
        }
    }
private:
    string name;
    bool hadValue;
    string previousValue;
};

struct PythonForceEvaluation {
    double energy;
    vector<Vec3> forces;
    bool callbackOnCallerThread;
};

void testForce() {
    class Computation : public PythonForceComputation {
        void compute(const State& state, double& energy, void* forces, bool forcesAreDouble) const {
            ASSERT_EQUAL(5.0, state.getParameters().at("a"));
            ASSERT_EQUAL(10.0, state.getParameters().at("b"));
            Vec3 a, b, c;
            state.getPeriodicBoxVectors(a, b, c);
            ASSERT_EQUAL(Vec3(2, 0, 0), a);
            ASSERT_EQUAL(Vec3(0.1, 2, 0), b);
            ASSERT_EQUAL(Vec3(0.1, 0.1, 2), c);
            energy = 25.0;
            int numParticles = state.getPositions().size();
            for (int i = 0; i < numParticles; i++) {
                Vec3 f = state.getPositions()[i]*2;
                if (forcesAreDouble)
                    ((Vec3*) forces)[i] = f;
                else {
                    ((float*) forces)[3*i] = (float) f[0];
                    ((float*) forces)[3*i+1] = (float) f[1];
                    ((float*) forces)[3*i+2] = (float) f[2];
                }
            }
        }
    };
    int numParticles = 5;
    System system;
    Vec3 a(2, 0, 0);
    Vec3 b(0.1, 2, 0);
    Vec3 c(0.1, 0.1, 2);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int i = 0; i < numParticles; i++) {
        system.addParticle(1.0);
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));
    }
    map<string, double> params;
    params["a"] = 5.0;
    params["b"] = 10.0;
    PythonForce* force = new PythonForce(new Computation(), params);
    ASSERT(!force->usesPeriodicBoundaryConditions());
    force->setUsesPeriodicBoundaryConditions(true);
    ASSERT(force->usesPeriodicBoundaryConditions());
    system.addForce(force);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    ASSERT_EQUAL_TOL(25.0, state.getPotentialEnergy(), 1e-6);
    for (int i = 0; i < numParticles; i++)
        ASSERT_EQUAL_VEC(2*positions[i], state.getForces()[i], 1e-6)

    // Check that force groups are handled correctly.

    ASSERT_EQUAL_TOL(25.0, context.getState(State::Energy, false, 1).getPotentialEnergy(), 1e-6);
    ASSERT_EQUAL_TOL(0.0, context.getState(State::Energy, false, 2).getPotentialEnergy(), 1e-6);
}

PythonForceEvaluation evaluatePythonForceExecutionMode(const vector<Vec3>& positions, const string& syncFlag) {
    class Computation : public PythonForceComputation {
    public:
        explicit Computation(thread::id callerThread) : callerThread(callerThread), callbackOnCallerThread(false) {
        }
        void compute(const State& state, double& energy, void* forces, bool forcesAreDouble) const override {
            callbackOnCallerThread = (this_thread::get_id() == callerThread);
            energy = 25.0;
            int numParticles = state.getPositions().size();
            for (int i = 0; i < numParticles; i++) {
                Vec3 f = state.getPositions()[i]*2;
                if (forcesAreDouble)
                    ((Vec3*) forces)[i] = f;
                else {
                    ((float*) forces)[3*i] = (float) f[0];
                    ((float*) forces)[3*i+1] = (float) f[1];
                    ((float*) forces)[3*i+2] = (float) f[2];
                }
            }
        }
        mutable bool callbackOnCallerThread;
    private:
        thread::id callerThread;
    };

    System system;
    for (int i = 0; i < (int) positions.size(); i++)
        system.addParticle(1.0);
    Computation* computation = new Computation(this_thread::get_id());
    system.addForce(new PythonForce(computation, map<string, double>()));
    ScopedEnvironmentVariable syncMode("OPENMM_PYTHONFORCE_SYNC", syncFlag);
    VerletIntegrator integrator(0.01);
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Energy | State::Forces);
    PythonForceEvaluation result;
    result.energy = state.getPotentialEnergy();
    result.forces = state.getForces();
    result.callbackOnCallerThread = computation->callbackOnCallerThread;
    return result;
}

void testSynchronousExecutionMode() {
    vector<Vec3> positions;
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(1, sfmt);
    for (int i = 0; i < 5; i++)
        positions.push_back(Vec3(genrand_real2(sfmt), genrand_real2(sfmt), genrand_real2(sfmt)));

    PythonForceEvaluation asyncMode = evaluatePythonForceExecutionMode(positions, "0");
    PythonForceEvaluation falseTextMode = evaluatePythonForceExecutionMode(positions, "FALSE");
    PythonForceEvaluation syncMode = evaluatePythonForceExecutionMode(positions, "1");

    ASSERT_EQUAL_TOL(25.0, asyncMode.energy, 1e-6);
    ASSERT_EQUAL_TOL(25.0, falseTextMode.energy, 1e-6);
    ASSERT_EQUAL_TOL(25.0, syncMode.energy, 1e-6);
    for (int i = 0; i < (int) positions.size(); i++) {
        ASSERT_EQUAL_VEC(2*positions[i], asyncMode.forces[i], 1e-6)
        ASSERT_EQUAL_VEC(2*positions[i], falseTextMode.forces[i], 1e-6)
        ASSERT_EQUAL_VEC(2*positions[i], syncMode.forces[i], 1e-6)
    }

    if (platform.getName() == "Reference") {
        ASSERT(asyncMode.callbackOnCallerThread);
        ASSERT(falseTextMode.callbackOnCallerThread);
        ASSERT(syncMode.callbackOnCallerThread);
    }
    else if (syncMode.callbackOnCallerThread) {
        ASSERT(!asyncMode.callbackOnCallerThread);
        ASSERT(!falseTextMode.callbackOnCallerThread);
    }
}

void runPlatformTests();

int main(int argc, char* argv[]) {
    try {
        initializeTests(argc, argv);
        testForce();
        testSynchronousExecutionMode();
        runPlatformTests();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
