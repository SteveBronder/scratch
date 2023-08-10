#include <iostream>
#include <vector>
#include <sundials/sundials_types.h>
#include <nvector/nvector_serial.h>
#include <arkode/oscodezs_arkstep.h>
#include <boost/math/special_functions/airy.hpp>

// Define the system of ODEs:
int airy_system(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    // Retrieve values from y vector:
    double u1 = NV_Ith_S(y, 0);
    double u2 = NV_Ith_S(y, 1);

    // Populate ydot:
    NV_Ith_S(ydot, 0) = u2;
    NV_Ith_S(ydot, 1) = t * u1;

    return 0;
}

int main() {
    // Define time span and steps:
    double T0 = 0.0;
    double Tf = 2.0;
    double dt = 0.01;

    // Initial conditions:
    double y0_data[2] = { boost::math::airy_ai(0.0), boost::math::airy_ai_prime(0.0) };

    // Create serial vector for storing solution:
    N_Vector y0 = N_VMake_Serial(2, y0_data);

    // Initialize ARKODE:
    void *arkode_mem = OSCStepCreate(airy_system, NULL, T0, y0);
    if (arkode_mem == NULL) {
        std::cerr << "Error initializing ARKODE" << std::endl;
        return -1;
    }

    // Specify tolerances:
    double reltol = 1.0e-6, abstol = 1.0e-8;
    OSCStepSStolerances(arkode_mem, reltol, abstol);

    // Integrate over time steps:
    double t = T0;
    while (t < Tf) {
        int flag = OSCStepEvolve(arkode_mem, t+dt, y0, &t, ARK_NORMAL);
        if (flag < 0) {
            std::cerr << "ARKStep failed at t=" << t << std::endl;
            return -1;
        }

        std::cout << "At t = " << t << ", y = [" << NV_Ith_S(y0, 0) << ", " << NV_Ith_S(y0, 1) << "]" << std::endl;
    }

    // Free memory:
    N_VDestroy(y0);
    OSCStepFree(&arkode_mem);

    return 0;
}
