#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "unused.hxx"
#include "output.hxx"

#include "reaction.hxx"
#include "radiation.hxx"

class ReactionIonchargedH2 : public Reaction {
public:
  ReactionIonchargedH2(Options *options) {
    TRACE("ReactionIonchargedH2(Options*)");

    // Ionisation energy cost
    OPTION(options, Ediss, 0.5);   /// bind energy between two atoms in a hydrogen molecule.

    bool diagnose;
    OPTION(options, diagnose, false);
    if (diagnose) {
      SAVE_REPEAT4(Riz, Eiz, Fiz, Siz);
    }
  }

  void updateSpecies(const SpeciesMap &species, BoutReal Tnorm,
                     BoutReal Nnorm, BoutReal UNUSED(Cs0), BoutReal Omega_ci) {
    // Get the species
    // Extract required variables
    auto &charged_M = *species.at("h2+");
    auto &electrons = *species.at("e");
//    auto &atoms     = *species.at('h')
    // Extract required variables
    Field3D Ne{electrons.N}, Te{electrons.T};
    Field3D Nm_c{charged_M.N}, Tm_c{charged_M.T}, Vm_c{charged_M.V};
        
    //Field3D Nm_c, Tm_c, Vm_c;
    //try {
    //  auto &charged_molecules = *species.at("h2+");
    //  Nm_c = charged_molecules.N;
    //  Tm_c = charged_molecules.T;
    //  Vm_c = charged_molecules.V;
    //} catch (const std::out_of_range &e) 
    //  throw BoutException("No 'h2+' species");
    //}


 
    Coordinates *coord = bout::globals::mesh->getCoordinates();

    Riz = 0.0;
    Eiz = 0.0;
    Fiz = 0.0;
    Siz = 0.0;
    
    CELL_AVERAGE(i,                        // Index variable
                 Riz.getRegion(RGN_NOBNDRY),  // Index and region (input)
                 coord,                    // Coordinate system (input)
                 weight,                   // Quadrature weight variable
                 Ne, Nm_c, Te, Tm_c, Vm_c) {     // Field variables

      BoutReal R =
        Ne * Nm_c * hydrogen.Ion_charged_H2(Te*Tnorm) * (Nnorm / Omega_ci);

      // Electron energy loss per ionisation
      Riz[i] += weight * (Ediss / Tnorm) * R;

      // Energy from neutral atom temperature
      Eiz[i] -= weight * (3. / 2) * Tm_c * R;

      // Friction due to ionisation
      Fiz[i] -= weight * Vm_c * R;

      // Plasma sink due to ionisation (negative)
      Siz[i] -= weight * R;
    }
  }
  
  SourceMap densitySources() override {
    return {{"h2+", Siz},  // Siz < 0 => H2+ source 
            {"h+", -2*Siz}};   // Siz < 0 => H+ sink
  }
  SourceMap momentumSources() {
    return {{"h2+",Fiz},
            {"h+",-0.5*Fiz}};  
  }
  SourceMap energySources() {
    return {{"h2+", Eiz}, // Eiz < 0 => Ion energy source
            {"e", -Riz},  // Electron energy into ionisation
            {"h+", -0.5*Eiz}};  // 
  }
  
  std::string str() const { return "charged molecule dissociation"; }
  
private:
  UpdatedRadiatedPower hydrogen; // Atomic rates
  Field3D Riz, Eiz, Fiz, Siz;

  BoutReal Ediss;  // Ionisation energy loss [eV]
};

namespace {
RegisterInFactory<Reaction, ReactionIonchargedH2> register_iz("ion_m");
}
