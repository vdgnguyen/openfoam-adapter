#include "Pressure.H"

using namespace Foam;

preciceAdapter::FSI::Pressure::Pressure(
    const Foam::fvMesh& mesh,
    const std::string solverType)
: ForceBase(mesh, solverType)
{
   
}

void preciceAdapter::FSI::Pressure::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    // Density boundary field
    tmp<volScalarField> trho(rho());
    const volScalarField::Boundary& rhob =
        trho().boundaryField();

    // Pressure boundary field
    const auto& pb = mesh_.lookupObject<volScalarField>("p").boundaryField();
    
    // For every boundary patch of the interface
    for (const label patchID : patchIDs_)
    {
        // For every cell of the patch
        forAll(pb[patchID], i)
        {
            // value in each cell
            double P = pb[patchID][i];
            double rho = rhob[patchID][i];
            if (solverType_.compare("incompressible") == 0)
            {
                buffer[i * dim] = -1*P*rho;
                for (unsigned int d = 1; d < dim; ++d)
                    buffer[i * dim+d] = 0;
                    
            }
            else if (solverType_.compare("compressible") == 0)
            {
                buffer[i * dim] = -1*P;
                for (unsigned int d = 1; d < dim; ++d)
                    buffer[i * dim +d] = 0;
            }
            else
            {
                FatalErrorInFunction
                    << "Pressure calculation does only support "
                    << "compressible or incompressible solver type."
                    << exit(FatalError);
            }
        }
    }
}

void preciceAdapter::FSI::Pressure::read(double* buffer, const unsigned int dim)
{
    this->readFromBuffer(buffer);
}

bool preciceAdapter::FSI::Pressure::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FSI::Pressure::getDataName() const
{
    return "Pressure";
}

Foam::tmp<Foam::vectorField> preciceAdapter::FSI::Pressure::getFaceVectors(const unsigned int patchID) const
{
    // face normal vectors
    return mesh_.boundary()[patchID].nf();
}

preciceAdapter::FSI::Pressure::~Pressure()
{

}
