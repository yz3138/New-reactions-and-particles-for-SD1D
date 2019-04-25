// Hydrogen charge exchange

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "field_factory.hxx"
#include "globals.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "radiation.hxx"
#include "reaction.hxx"

using bout::globals::mesh;

class ReactionHydrogenDissociation : public Reaction {
public:
  ReactionHydrogenDissociation(Options *options) {
    AUTO_TRACE();

    OPTION(options, Ediss, 0.5);   // eV, The bind energy between two atoms in a hydrogen molecule

    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT4(Sn,En,Ee,Fn);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal UNUSED(Cs0), BoutReal Omega_ci) {

    // Get the species
    auto &molecules = *species.at("h2");
    auto &electrons = *species.at("e");
//    auto &atoms     = *species.at('h')
    // Extract required variables
    Field3D Ne{electrons.N}, Te{electrons.T};
    Field3D Nm{molecules.N}, Tm{molecules.T}, Vm{molecules.V};
//    Field3D Nn{atoms.N}, Tn{atoms.T}, Vn{atoms.V};
    Coordinates *coord = mesh->getCoordinates();

    Sn = 0.0;
    En = 0.0;
    Fn = 0.0;
    Ee = 0.0;
//    Sn_1=0.0;
    CELL_AVERAGE(i,                          // Index variable
                 Sn.getRegion(RGN_NOBNDRY), // Index and region (input)
                 //Sn_1.getRegion(RGN_NOBNDRY),
                 coord,                      // Coordinate system (input)
                 weight,                     // Quadrature weight variable
                 Ne, Nm, Te, Tm, Vm) {   // Field variables

      BoutReal R = Ne * Nm * hydrogen.dissociation_ion(Te*Tnorm) * (Nnorm / Omega_ci);  // dissociation ionization (DI)
//      BoutReal R1 = Nn * Nm * hydrogen.Diss_H_H2(Te*Tnorm) * (Nnorm / Omega_ci);
// Dissociation by collision between H and H2

                            // ~~~~~  Rate here in SI units m^3/s

      // Ecx is energy transferred from ions to neutrals
      Sn[i] += weight * R;  //Source term from DI and Diss_H_H2
//      Sn_1[i] += weight * R1;

      // Energy loss per dissociation
      Ee[i] += weight * (Ediss/Tnorm) * R;
      // Energy from molecule temperature
      En[i] -= weight * (3./2.) * Tm * R;
      Fn[i] -= weight * Vm * R;
   
    }
  }
  SourceMap densitySources() override {
    return {{"h", Sn},   ///????  why is not '2*Sn' for h, and 'Sn' for h2????
            {"h2", -Sn},
            {"h+",Sn}}; // Redistribute neutrals
  }
  SourceMap momentumSources() {
    return {{"h2",Fn},
            {"h",-0.5*Fn},
            {"h+",-0.5*Fn}};
  }
  SourceMap energySources() {
    return {{"h2",En},
            {"e",-Ee},
            {"h",-0.5*En},
            {"h+",-0.5*En}         };
  }

  std::string str() const { return "Hydrogen dissociation"; }

private:
  Field3D Sn,En,Ee,Fn; // Source of atoms
  BoutReal Ediss;
  UpdatedRadiatedPower hydrogen; // Atomic rates
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenDissociation> register_rc("hydrogen_diss");
}
