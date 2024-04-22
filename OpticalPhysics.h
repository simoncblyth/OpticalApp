#pragma once

#include "G4VUserPhysicsList.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpAbsorption.hh"


struct OpticalPhysics : public G4VUserPhysicsList
{
    const char* config ; 

    OpticalPhysics(const char* config); 

    void ConstructParticle(); 
    void ConstructProcess() ; 
    void ConstructOp();
};

OpticalPhysics::OpticalPhysics(const char* _config)
    :
    config(_config ? strdup(_config) : nullptr )
{
} 

inline void OpticalPhysics::ConstructParticle()
{
    G4OpticalPhoton::OpticalPhotonDefinition();
}
inline void OpticalPhysics::ConstructProcess()
{
    AddTransportation();
    ConstructOp();    
}

inline void OpticalPhysics::ConstructOp()
{

    G4VProcess* boundary = new G4OpBoundaryProcess();   
    G4VProcess* absorption = new G4OpAbsorption();   

    auto particleIterator=GetParticleIterator();
    particleIterator->reset();
    while( (*particleIterator)() )
    {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (particleName == "opticalphoton")
        {
            std::stringstream ss;  
            ss.str(config)  ;
            std::string cfg;
            while (std::getline(ss, cfg, ',')) 
            {   
                if(cfg.empty()) continue ;   
                const char* _cfg = cfg.c_str() ;  
                if(      strcmp(_cfg,"G4OpAbsorption")==0) pmanager->AddDiscreteProcess(absorption);
                else if( strcmp(_cfg,"G4OpBoundaryProcess")==0) pmanager->AddDiscreteProcess(boundary);
                else
                {
                   std::cerr << "OpticalPhysics::ConstructOp FATAL unhandled [" << ( _cfg ? _cfg : "-" ) << "]\n" ;
                }
            }

        }
    }
}

