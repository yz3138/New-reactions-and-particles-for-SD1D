#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "globals.hxx"
#include "unused.hxx"
#include "output.hxx"

#include "reaction.hxx"
#include "radiation.hxx"

class ReactionH2Ionisation : public Reaction {
public:
  ReactionH2Ionisation(Options *options) {
    TRACE("ReactionH2Ionisation(Options*)");

    // Ionisation energy cost
    OPTION(options, Eionize, 30);   /// may different for H2->H2+ ???

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
    auto &molecules = *species.at("h2");
    auto &electrons = *species.at("e");
//    auto &atoms     = *species.at('h')
    // Extract required variables
    Field3D Ne{electrons.N}, Te{electrons.T};
    Field3D Nm{molecules.N}, Tm{molecules.T}, Vm{molecules.V};
        
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
                 Ne, Nm, Te, Tm, Vm) {     // Field variables

      BoutReal R =
        Ne * Nm * hydrogen.Ion_H2(Te*Tnorm) * (Nnorm / Omega_ci);

      // Electron energy loss per ionisation
      Riz[i] += weight * (Eionize / Tnorm) * R;

      // Energy from neutral atom temperature
      Eiz[i] -= weight * (3. / 2) * Tm * R;

      // Friction due to ionisation
      Fiz[i] -= weight * Vm * R;

      // Plasma sink due to ionisation (negative)
      Siz[i] -= weight * R;
    }
  }
  
  SourceMap densitySources() override {
    return {{"h2+", -Siz},  // Siz < 0 => H2+ source 
            {"h2", Siz}};   // Siz < 0 => H2 sink
  }
  SourceMap momentumSources() {
    return {{"h2+", -Fiz}, // charged molecule momentum source
            {"h2", Fiz}};  // Neutral molecule momentum sink
  }
  SourceMap energySources() {
    return {{"h2+", -Eiz}, // Eiz < 0 => Ion energy source
            {"e", -Riz},  // Electron energy into ionisation
            {"h2", Eiz}};  // Neutral atom energy sink
  }
  
  std::string str() const { return "molecule ionisation"; }
  
private:
  UpdatedRadiatedPower hydrogen; // Atomic rates
  Field3D Riz, Eiz, Fiz, Siz;

  BoutReal Eionize;  // Ionisation energy loss [eV]
};

namespace {
RegisterInFactory<Reaction, ReactionH2Ionisation> register_iz("ion_h2");
}
