#ifndef FSI_PRESSURE_H
#define FSI_PRESSURE_H

#include "ForceBase.H"

namespace preciceAdapter
{
namespace FSI
{

//- Class that writes and reads stress [N/m^2]:
// This is essentially a force (in spatial coordinates) scaled by the deformed
// cell face. Thus, a consistent quantity. Calculation concept has been copied
// from the force module, but the scaled version here is commonly used in FEM
// applications.
class Pressure : public ForceBase
{

public:
    //- Constructor
    Pressure(
        const Foam::fvMesh& mesh,
        const std::string solverType);

    //- Write the stress values into the buffer
    void write(double* buffer, bool meshConnectivity, const unsigned int dim) final;

    //- Read the stress values from the buffer
    void read(double* buffer, const unsigned int dim) final;

    bool isLocationTypeSupported(const bool meshConnectivity) const final;

    //- Get the name of the current data field
    std::string getDataName() const final;

    //- not used
    Foam::tmp<Foam::vectorField> getFaceVectors(const unsigned int patchID) const final;

    //- Destructor
    ~Pressure();
};
}
}

#endif
