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

class ReactionHydrogenDissH : public Reaction {
public:
  ReactionHydrogenDissH(Options *options) {
    AUTO_TRACE();


    OPTION(options, Ediss, 0.5);   // eV, The bind energy between two atoms in a hydrogen molecule


    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT(Sn,En,Em,Fn);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm, BoutReal Nnorm,
                     BoutReal UNUSED(Cs0), BoutReal Omega_ci) {

    // Get the species
    auto &atoms     = *species.at("h");
    auto &molecules = *species.at("h2");
//    auto &electrons = *species.at("e");
//    auto &atoms     = *species.at('h')
    // Extract required variables
//    Field3D Ne{electrons.N}, Te{electrons.T};
    Field3D Nm{molecules.N}, Tm{molecules.T}, Vm{molecules.V};
    Field3D Nn{atoms.N}, Tn{atoms.T}, Vn{atoms.V};
    Coordinates *coord = mesh->getCoordinates();

    Sn = 0.0;
    En = 0.0;
    Fn = 0.0;
    Em = 0.0;
//    Sn_1=0.0;
    CELL_AVERAGE(i,                          // Index variable
                 Sn.getRegion(RGN_NOBNDRY), // Index and region (input)
                 //Sn_1.getRegion(RGN_NOBNDRY),
                 coord,                      // Coordinate system (input)
                 weight,                     // Quadrature weight variable
                 Nm, Nn, Tm, Tn, Vm, Vn) {   // Field variables

//      BoutReal R = Ne * Nm * hydrogen.dissociation_ion(Te * Tnorm) * (Nnorm / Omega_ci);  // dissociation ionization (DI)
      BoutReal R = Nn * Nm * hydrogen.Diss_H_H2(Tn*Tnorm) * (Nnorm / Omega_ci);
// Dissociation by collision between H and H2

                            // ~~~~~  Rate here in SI units m^3/s

      // Ecx is energy transferred from ions to neutrals
      Sn[i] += weight * R;  //Source term from DI and Diss_H_H2
//      Sn_1[i] += weight * R1;
      En[i] -= weight * (Ediss/Tnorm) * R;
      Em[i] -= weight * (3./2.) * Tm * R; // energy transferred from molecule to atom
      Fn[i] -= weight * (Vn-Vm) * R;
 
    }
  }
  SourceMap densitySources() override {
    return {{"h", Sn},   ///????  why is not '2*Sn' for h, and 'Sn' for h2????
            {"h2", -0.5*Sn}}; // Redistribute neutrals
  }
  SourceMap momentumSources() {
    return {{"h", Fn},
            {"h2",-Fn}};
  }
  SourceMap energySources() {
    return {{"h", En-0.5*Em},   ///  H2->H+H, so the energy also splitted by 2
            {"h2",Em}};
  }

  std::string str() const { return "Hydrogen dissociation_DissH"; }

private:
  Field3D Sn,En,Fn,Em; // Source of atoms
  UpdatedRadiatedPower hydrogen; // Atomic rates
  BoutReal Ediss;
};

namespace {
RegisterInFactory<Reaction, ReactionHydrogenDissH> register_rc("hydrogen_dissh");
} /// cannot register with upper cases!!!!!

