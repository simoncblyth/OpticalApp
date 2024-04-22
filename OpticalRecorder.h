#pragma once

#include <algorithm>
#include <map>
#include "G4OpticalPhoton.hh"
#include "G4OpBoundaryProcess.hh"
#include "np.h"

#define FFS(x)   (ffs(x))

// OpticksPhoton.h 
enum
{
    CERENKOV          = 0x1 <<  0,    
    SCINTILLATION     = 0x1 <<  1,    
    MISS              = 0x1 <<  2,  
    BULK_ABSORB       = 0x1 <<  3,  
    BULK_REEMIT       = 0x1 <<  4,  
    BULK_SCATTER      = 0x1 <<  5,  
    SURFACE_DETECT    = 0x1 <<  6,  
    SURFACE_ABSORB    = 0x1 <<  7,  
    SURFACE_DREFLECT  = 0x1 <<  8,  
    SURFACE_SREFLECT  = 0x1 <<  9,  
    BOUNDARY_REFLECT  = 0x1 << 10, 
    BOUNDARY_TRANSMIT = 0x1 << 11, 
    TORCH             = 0x1 << 12, 
    NAN_ABORT         = 0x1 << 13, 
    EFFICIENCY_CULL    = 0x1 << 14, 
    EFFICIENCY_COLLECT = 0x1 << 15, 
    __NATURAL         = 0x1 << 16, 
    __MACHINERY       = 0x1 << 17, 
    __EMITSOURCE      = 0x1 << 18, 
    PRIMARYSOURCE     = 0x1 << 19, 
    GENSTEPSOURCE     = 0x1 << 20, 
    DEFER_FSTRACKINFO = 0x1 << 21
}; 


union uif64_t {
        uint64_t  u ; 
        int64_t   i ; 
        double    f ; 
};  


struct OpticalRecorder
{
    static constexpr const unsigned BITS = 4 ;
    static constexpr const uint64_t MASK = ( 0x1ull << BITS ) - 1ull ;

    static constexpr const char* fGeomBoundary_  = "fGeomBoundary" ;
    static constexpr const char* fUndefined_ = "fUndefined" ;
    static constexpr const char* fERROR_ = "fERROR" ;
    static const char* StepStatus(unsigned status) ; 

    static constexpr const char* Undefined_ = "Undefined" ;
    static constexpr const char* Transmission_ = "Transmission" ;
    static constexpr const char* FresnelRefraction_ = "FresnelRefraction" ;
    static constexpr const char* FresnelReflection_  = "FresnelReflection" ;
    static constexpr const char* TotalInternalReflection_ = "TotalInternalReflection" ;
    static constexpr const char* LambertianReflection_ = "LambertianReflection" ;
    static constexpr const char* LobeReflection_ = "LobeReflection" ;
    static constexpr const char* SpikeReflection_ = "SpikeReflection" ;
    static constexpr const char* BackScattering_ = "BackScattering" ;
    static constexpr const char* Absorption_ = "Absorption" ;
    static constexpr const char* Detection_ = "Detection" ;
    static constexpr const char* NotAtBoundary_ = "NotAtBoundary" ;
    static constexpr const char* SameMaterial_ = "SameMaterial" ;
    static constexpr const char* StepTooSmall_ = "StepTooSmall" ;
    static constexpr const char* NoRINDEX_ = "NoRINDEX" ;
    static constexpr const char* Other_ = "Other" ;   
    static const char* OpBoundaryProcessStatus(unsigned status); 


    static constexpr const char* _ZERO              = "  " ;
    static constexpr const char* _CERENKOV          = "CK" ;
    static constexpr const char* _SCINTILLATION     = "SI" ; 
    static constexpr const char* _TORCH             = "TO" ; 
    static constexpr const char* _MISS              = "MI" ;
    static constexpr const char* _BULK_ABSORB       = "AB" ;
    static constexpr const char* _BULK_REEMIT       = "RE" ;
    static constexpr const char* _BULK_SCATTER      = "SC" ;
    static constexpr const char* _SURFACE_DETECT    = "SD" ;
    static constexpr const char* _SURFACE_ABSORB    = "SA" ;
    static constexpr const char* _SURFACE_DREFLECT  = "DR" ;
    static constexpr const char* _SURFACE_SREFLECT  = "SR" ;
    static constexpr const char* _BOUNDARY_REFLECT  = "BR" ;
    static constexpr const char* _BOUNDARY_TRANSMIT = "BT" ;
    static constexpr const char* _NAN_ABORT         = "NA" ;
    static constexpr const char* _EFFICIENCY_COLLECT = "EC" ;
    static constexpr const char* _EFFICIENCY_CULL    = "EX" ;
    static constexpr const char* _BAD_FLAG           = "XX" ;
    static const char* Flag(unsigned flag); 
    static std::string FlagSequence(uint64_t seqhis ); 
    static const char* FlagElement(uint64_t seqhis, unsigned j);


    template <typename T>
    static T* GetOpBoundaryProcess();

    void writePoint( const G4StepPoint* point, unsigned flag ); 
    uint64_t getSeq(int _trk_idx) const ; 
 
    static void WritePoint( double* p , const G4StepPoint* point, unsigned flag ); 
    static std::string Desc( const double* p, int num ); 
 
    void BeginOfRunAction(const G4Run* run);
    void EndOfRunAction(const G4Run* run);

    void BeginOfEventAction(const G4Event* evt);
    void EndOfEventAction(const G4Event* evt);

    void PreUserTrackingAction(const G4Track* trk);
    void PostUserTrackingAction(const G4Track* trk);

    static std::string DescSeq( const std::map<uint64_t,uint64_t>& _seqmap ); 
    static void GetSeqName( std::vector<std::string>& seqname, const std::map<uint64_t,uint64_t>& _seqmap ); 
    std::string descSeq() const ; 

    static bool StartsWith( const char* s, const char* q); 
    static uint64_t FindSeq(const std::map<uint64_t,uint64_t>& _seqmap, const char* seq ); 
    static uint64_t CountNibbles(uint64_t x); 

    std::string descSpeed(const char* seq) const ; 
    std::string descSpeedAll() const ; 

    static int TrackIdx( const G4Track* trk ); 

    static unsigned PointFlag( const G4StepPoint* point ); 
    static unsigned BoundaryFlag(unsigned status);  // BT BR NA SA SD SR DR

    void UserSteppingAction(const G4Step* step);


    static bool Valid(int trk_idx, int point_idx);
    static const double* GetRecord(const double* _pp, int _trk_idx, int _point_idx); 

    const double* getRecord(int _point_idx) const ;
    void recordPoint( const G4StepPoint* point ); 

    static unsigned PointFlag(const double* a); 
    static double Wavelength( const double* a ); 
    static double DeltaTime( const double* a, const double* b );
    static double DeltaPos( const double* a, const double* b );
    static double GetSpeed(const double* _pp, int _trk_idx, int _point_idx ); 
    static void  GetPointSpeed( std::vector<double>& speeds, const double* _pp, uint64_t seqhis, int _point_idx ); 

    double getSpeed(int _point_idx) const ; 
    std::string descPoint(int _point_idx) const ; 


    static uint64_t GetHistory(const double* _pp,  int _trk_idx); 
    static std::string DescHistory(const double* _pp, int _trk_idx) ; 
    std::string descHistory(int _trk_idx) const ; 
    uint64_t    getHistory(int _trk_idx) const ; 

    static constexpr const int MAX_PHOTON = 100000 ; 
    static constexpr const int MAX_POINT  = 10 ; 
    static constexpr const int MAX_VALUE  = MAX_PHOTON*MAX_POINT*16  ; 
  
    OpticalRecorder(); 
    void clear(); 
    void alloc(); 

    int trk_idx = 0  ; 
    int point_idx = 0  ; 

    double* pp = nullptr ; 
    uint64_t* qq = nullptr ; 
    std::string desc ; 

    std::map<uint64_t,uint64_t> seqmap ; 

}; 

inline OpticalRecorder::OpticalRecorder(){}

inline void OpticalRecorder::clear()
{
    delete [] pp ; 
    pp = nullptr ; 

    delete [] qq ; 
    qq = nullptr ; 
}
inline void OpticalRecorder::alloc()
{
    pp = new double[MAX_VALUE] ; 
    for(int i=0 ; i < MAX_VALUE ; i++) pp[i] = 0. ; 

    qq = new uint64_t[MAX_PHOTON*2*2] ; 
    for(int i=0 ; i < MAX_PHOTON*2*2 ; i++) qq[i] = 0ull ; 
}


// U4StepStatus::Name
inline const char* OpticalRecorder::StepStatus(unsigned status)
{
    const char* str = nullptr ; 
    switch(status)
    {   
        case fGeomBoundary:           str=fGeomBoundary_           ;break; 
        case fUndefined:              str=fUndefined_              ;break; 
        default:                      str=fERROR_                  ;break;
    }   
    return str ; 
}

// U4OpBoundaryProcessStatus::Name
inline const char* OpticalRecorder::OpBoundaryProcessStatus(unsigned status)
{
    const char* str = nullptr ; 
    switch(status)
    {   
       case Undefined:                str = Undefined_ ; break ; 
       case Transmission:             str = Transmission_ ; break ; 
       case FresnelRefraction:        str = FresnelRefraction_ ; break ; 
       case FresnelReflection:        str = FresnelReflection_ ; break ;
       case TotalInternalReflection:  str = TotalInternalReflection_ ; break ;
       case LambertianReflection:     str = LambertianReflection_ ; break ;
       case LobeReflection:           str = LobeReflection_ ; break ;
       case SpikeReflection:          str = SpikeReflection_ ; break ;
       case BackScattering:           str = BackScattering_ ; break ;
       case Absorption:               str = Absorption_ ; break ;
       case Detection:                str = Detection_ ; break ;
       case NotAtBoundary:            str = NotAtBoundary_ ; break ;
       case SameMaterial:             str = SameMaterial_ ; break ;
       case StepTooSmall:             str = StepTooSmall_ ; break ;
       case NoRINDEX:                 str = NoRINDEX_ ; break ;
       default:                       str = Other_ ; break ;
    }
    return str ;
}

// OpticksPhoton::Flag
const char* OpticalRecorder::Flag(const unsigned int flag)
{
    const char* str = 0 ;
    switch(flag)
    {
        case 0:                str=_ZERO;break;
        case CERENKOV:         str=_CERENKOV ;break;
        case SCINTILLATION:    str=_SCINTILLATION ;break;
        case MISS:             str=_MISS ;break;
        case BULK_ABSORB:      str=_BULK_ABSORB ;break;
        case BULK_REEMIT:      str=_BULK_REEMIT ;break;
        case BULK_SCATTER:     str=_BULK_SCATTER ;break;
        case SURFACE_DETECT:   str=_SURFACE_DETECT ;break;
        case SURFACE_ABSORB:   str=_SURFACE_ABSORB ;break;
        case SURFACE_DREFLECT: str=_SURFACE_DREFLECT ;break;
        case SURFACE_SREFLECT: str=_SURFACE_SREFLECT ;break;
        case BOUNDARY_REFLECT: str=_BOUNDARY_REFLECT ;break;
        case BOUNDARY_TRANSMIT:str=_BOUNDARY_TRANSMIT ;break;
        case TORCH:            str=_TORCH ;break;
        case NAN_ABORT:        str=_NAN_ABORT ;break;
        case EFFICIENCY_CULL:    str=_EFFICIENCY_CULL ;break;
        case EFFICIENCY_COLLECT: str=_EFFICIENCY_COLLECT ;break;
        default:                 str=_BAD_FLAG  ;
    }
    return str;
}

// OpticksPhoton::FlagSequence
std::string OpticalRecorder::FlagSequence(uint64_t seqhis )
{
    std::stringstream ss ;
    for(unsigned j=0 ; j < 16   ; j++) ss << FlagElement(seqhis,j) << " " ; 
    return ss.str();
}

const char* OpticalRecorder::FlagElement(uint64_t seqhis, unsigned j)
{
    uint64_t f = (seqhis >> j*4) & 0xF ; 
    unsigned flg = f == 0 ? 0 : 0x1 << (f - 1) ; 
    return Flag(flg) ; 
}


// U4OpBoundaryProcess::Get
template <typename T>
inline T* OpticalRecorder::GetOpBoundaryProcess()
{           
    T* bp = nullptr ;
    G4ProcessManager* mgr = G4OpticalPhoton::OpticalPhoton()->GetProcessManager() ;
    assert(mgr); 

    G4int pmax = mgr ? mgr->GetPostStepProcessVector()->entries() : 0 ;
    G4ProcessVector* pvec = mgr ? mgr->GetPostStepProcessVector(typeDoIt) : nullptr ;

    for (int i=0; i < pmax ; i++)
    {
        G4VProcess* p = (*pvec)[i];
        T* t = dynamic_cast<T*>(p);
        if(t)
        {
            bp = t ;
            break;
        }
    }
    return bp ;
}

// U4StepPoint::Update
void OpticalRecorder::writePoint( const G4StepPoint* point, unsigned flag )
{
    const double* p = getRecord(point_idx); 
    if(!p) return ; 

    WritePoint( const_cast<double*>(p), point, flag ); 

    uint64_t& q0 = qq[4*trk_idx+0] ; 

    // sseq::add_nibble
    unsigned shift = 4*point_idx ; 
    q0 |= (( FFS(flag) & MASK ) << shift );  
    
}

uint64_t OpticalRecorder::getSeq(int _trk_idx) const
{
    return qq[4*_trk_idx+0] ; 
}


void OpticalRecorder::WritePoint( double* p, const G4StepPoint* point, unsigned flag )
{
    const G4ThreeVector& pos = point->GetPosition();
    const G4ThreeVector& mom = point->GetMomentumDirection();
    const G4ThreeVector& pol = point->GetPolarization();

    G4double time = point->GetGlobalTime();
    G4double energy = point->GetKineticEnergy();
    G4double wavelength = CLHEP::h_Planck*CLHEP::c_light/energy ;
    
    p[0] = pos.x(); 
    p[1] = pos.y(); 
    p[2] = pos.z(); 
    p[3] = time/ns ; 

    p[4] = mom.x(); 
    p[5] = mom.y(); 
    p[6] = mom.z(); 
    p[7] = 0. ; 

    p[8] = pol.x();
    p[9] = pol.y();
    p[10] = pol.z();
    p[11] = wavelength/nm ;


    uif64_t uif ; 
    uif.u = flag ;  

    p[12] = 0. ; 
    p[13] = 0. ; 
    p[14] = 0. ; 
    p[15] = uif.f ; 
}



unsigned OpticalRecorder::PointFlag(const double* a)
{
    if(a == nullptr) return 0 ; 
    uif64_t uif ; 
    uif.f = a[15] ;  
    return unsigned(uif.u) ; 
}

double OpticalRecorder::Wavelength( const double* a )
{
    return a ? a[11] : 0. ; 
} 

double OpticalRecorder::DeltaTime( const double* a, const double* b )
{
    return a && b ? b[3] - a[3] : -1 ; 
}
double OpticalRecorder::DeltaPos( const double* a, const double* b )
{
    if( a == nullptr || b == nullptr ) return -1 ; 
    G4ThreeVector dpos( b[0] - a[0], b[1] - a[1], b[2] - a[2] ); 
    return dpos.mag() ; 
}







std::string OpticalRecorder::Desc( const double* p, int num )
{
    std::stringstream ss ; 
    assert( num == 16 ); 
    for(int i=0 ; i < num ; i++) 
        ss  
            << ( i % 4 == 0 && num > 4 ? ".\n" : "" ) 
            << " " << std::fixed << std::setw(10) << std::setprecision(4) << p[i] 
            << ( i == num-1 && num > 4 ? ".\n" : "" ) 
            ;   

    std::string str = ss.str(); 
    return str ; 
}



void OpticalRecorder::BeginOfRunAction(const G4Run* ){         std::cout << "OpticalRecorder::BeginOfRunAction\n" ;    }
void OpticalRecorder::EndOfRunAction(const G4Run* ){           std::cout << "OpticalRecorder::EndOfRunAction\n" ;  }


void OpticalRecorder::BeginOfEventAction(const G4Event* evt)
{  
    std::cout << "OpticalRecorder::BeginOfEventAction evt " << evt->GetEventID() << "\n" ;
    alloc();  
} 

void OpticalRecorder::EndOfEventAction(const G4Event* evt)
{    
    int eid = evt->GetEventID() ; 
    std::cout << "OpticalRecorder::EndOfEventAction eid " << eid << "\n" ; 

    const char* FOLD = getenv("FOLD") ; 

    std::vector<int> pp_shape = { MAX_PHOTON, MAX_POINT, 4, 4 } ; 
    np::Write( FOLD, "record.npy",  pp_shape, pp, "<f8" ); 

    std::vector<int> qq_shape = { MAX_PHOTON, 2, 2 } ; 
    np::Write( FOLD, "seq.npy",  qq_shape, qq, "<u8" ); 

    if(!desc.empty()) np::WriteString( FOLD, "NPFold_meta.txt", desc.c_str() ); 

    std::cout << descSeq() ; 
    //std::cout << descSpeed("TO BT SA") ; 
    //std::cout << descSpeed("TO BR BT SA") ; 
    //std::cout << descSpeed("TO BR BR BT SA") ; 

    std::cout << descSpeedAll() ; 

    clear(); 
}
void OpticalRecorder::PreUserTrackingAction(const G4Track* trk )
{ 
    trk_idx = TrackIdx(trk) ;
    point_idx = 0 ;
 
    if(0) std::cout 
        << "OpticalRecorder::PreUserTrackingAction"
        << " trk_idx :" << trk_idx 
        << " point_idx:" << point_idx 
        << "\n" 
        ;  
}
void OpticalRecorder::PostUserTrackingAction(const G4Track* trk)
{
    assert( TrackIdx(trk) == trk_idx ); 
 
    uint64_t seq = getHistory(trk_idx); 
    if(0) std::cout 
        << "OpticalRecorder::PostUserTrackingAction"
        << " trk_idx :" << trk_idx 
        << " point_idx:" << point_idx 
        << " descHistory: " << descHistory(trk_idx)
        << " seq: " << std::hex << seq << std::dec
        << "\n" 
        ; 


   if(seqmap.count(seq)==0) seqmap[seq] = 1 ; 
   else                     seqmap[seq] += 1 ;  
}


std::string OpticalRecorder::DescSeq( const std::map<uint64_t,uint64_t>& _seqmap )
{
    std::stringstream ss ;
    ss << "[OpticalRecorder::DescSeq\n" ;  
    typedef std::map<uint64_t,uint64_t> MUU ; 
    for(MUU::const_iterator it=_seqmap.begin() ; it!=_seqmap.end() ; it++) ss 
        << std::setw(10) << it->second
        << " : "
        << FlagSequence(it->first )
        << "\n"
        ;
    ss << "]OpticalRecorder::DescSeq\n" ;  
    return ss.str(); 
}

void OpticalRecorder::GetSeqName( std::vector<std::string>& seqname, const std::map<uint64_t,uint64_t>& _seqmap ) 
{
    typedef std::map<uint64_t,uint64_t> MUU ; 
    for(MUU::const_iterator it=_seqmap.begin() ; it!=_seqmap.end() ; it++) 
    {
        std::string nm = FlagSequence(it->first ) ;
        seqname.push_back(nm);    
    }
} 

std::string OpticalRecorder::descSeq() const 
{
    return DescSeq(seqmap); 
}


// sstr::StartsWith
bool OpticalRecorder::StartsWith( const char* s, const char* q)   // static
{
    return s && q && strlen(q) <= strlen(s) && strncmp(s, q, strlen(q)) == 0 ; 
}

uint64_t OpticalRecorder::FindSeq(const std::map<uint64_t,uint64_t>& _seqmap, const char* seq ) // static
{
    uint64_t seqhis = 0 ; 
    typedef std::map<uint64_t,uint64_t> MUU ;
    for(MUU::const_iterator it=_seqmap.begin() ; it!=_seqmap.end() ; it++) 
    {
        std::string nm = FlagSequence(it->first ) ;
        if(StartsWith(nm.c_str(), seq)) seqhis = it->first ;               
    }
    return seqhis ; 
}

// SBit::count_nibbles
uint64_t OpticalRecorder::CountNibbles(uint64_t x)
{
    x |= x >> 1 ; 
    x |= x >> 2 ; 
    x &= 0x1111111111111111ull ; 

    x = (x + (x >> 4)) & 0xF0F0F0F0F0F0F0Full ; 

    uint64_t count = (x * 0x101010101010101ull) >> 56 ; 
    return count ; 
}


std::string OpticalRecorder::descSpeed(const char* seq) const 
{
    uint64_t seqhis = FindSeq(seqmap, seq ); 
    uint64_t nibs = CountNibbles(seqhis); 
    std::stringstream ss ;
    ss << "[OpticalRecorder::descSpeed\n" ;  
    ss << " seq[" << ( seq ? seq : "-" ) << "]" ; 
    ss << " seqhis[" << std::hex << seqhis << std::dec << "]" ;  
    ss << " nibs[" << nibs << "]\n" ; 
    for(int i=1 ; i < int(nibs) ; i++)
    {
        int _point_idx = i ; 
        std::vector<double> speeds ; 
        GetPointSpeed( speeds, pp, seqhis, _point_idx ); 
        auto minmax = std::minmax_element( speeds.begin(), speeds.end() ); 
        ss
           << FlagElement(seqhis, i-1 )
           << " -> "
           << FlagElement(seqhis, i )
           << " speeds.len/min/max[" 
           << speeds.size() 
           << "/"  
           << *minmax.first 
           << "/"  
           << *minmax.second 
           << "]\n" 
           ;
    }
    ss << "]OpticalRecorder::descSpeed\n" ;  
    return ss.str(); 
}

std::string OpticalRecorder::descSpeedAll() const 
{
    std::stringstream ss ;
    typedef std::map<uint64_t,uint64_t> MUU ;
    for(MUU::const_iterator it=seqmap.begin() ; it!=seqmap.end() ; it++) 
    {
        std::string nm = FlagSequence(it->first ) ;
        ss << descSpeed(nm.c_str()) ; 
    }
    return ss.str(); 
}


int OpticalRecorder::TrackIdx( const G4Track* track )
{
   return track->GetTrackID() - 1 ;  // 0-based unlike 1-based TrackID  
}

// U4StepPoint::Flag
unsigned OpticalRecorder::PointFlag( const G4StepPoint* point )
{
    G4StepStatus status = point->GetStepStatus()  ;

    // U4StepPoint::ProcessDefinedStepType
    const G4VProcess* proc = point->GetProcessDefinedStep() ;
    const char* procName = proc ? proc->GetProcessName().c_str() : nullptr  ; 

    unsigned flag = 0 ; 

    if( status == fPostStepDoItProc && strcmp(procName, "OpAbsorption") == 0 )
    {
        flag = BULK_ABSORB ;
    }
    else if( status == fPostStepDoItProc && strcmp(procName, "OpRayleigh") == 0 )
    {
        flag = BULK_SCATTER ;
    }    
    else if( status == fGeomBoundary && strcmp(procName, "Transportation") == 0  )
    {
        G4OpBoundaryProcess* bp = GetOpBoundaryProcess<G4OpBoundaryProcess>(); 
        G4OpBoundaryProcessStatus bp_status = bp ? bp->GetStatus() : Undefined ;
        flag = BoundaryFlag( bp_status ); 
    }
    else if( status == fWorldBoundary && strcmp(procName, "Transportation") == 0  )
    {
        flag = MISS ;  
    }
    return flag  ; 
}

// U4StepPoint::BoundaryFlag
unsigned OpticalRecorder::BoundaryFlag(unsigned status) // BT BR NA SA SD SR DR 
{
    unsigned flag = 0 ; 
    switch(status)
    {   
        case FresnelRefraction:
        case SameMaterial:
        case Transmission:
                               flag=BOUNDARY_TRANSMIT;
                               break;
        case TotalInternalReflection:
        case       FresnelReflection:
                               flag=BOUNDARY_REFLECT;
                               break;
        case StepTooSmall:
                               flag=NAN_ABORT;
                               break;
        case Absorption:
                               flag=SURFACE_ABSORB ; 
                               break;
        case Detection:
                               flag=SURFACE_DETECT ; 
                               break;
        case SpikeReflection:
                               flag=SURFACE_SREFLECT ; 
                               break;
        case LobeReflection:
        case LambertianReflection:
                               flag=SURFACE_DREFLECT ; 
                               break;
        case NoRINDEX:
                               flag=SURFACE_ABSORB ;
                               //flag=NAN_ABORT;
                               break;
        default:
                               flag = 0 ; 
                               break;
    }
    return flag ; 
}


void OpticalRecorder::UserSteppingAction(const G4Step* step )
{    
    const G4Track* track = step->GetTrack();
    assert( trk_idx == TrackIdx(track) ); 

    const G4StepPoint* pre = step->GetPreStepPoint() ;
    const G4StepPoint* post = step->GetPostStepPoint() ;

    if(point_idx == 0) recordPoint(pre) ; 
    recordPoint(post); 
}

bool OpticalRecorder::Valid(int _trk_idx, int _point_idx)
{
    return 
          _trk_idx > -1 && _trk_idx < MAX_PHOTON 
          && 
          _point_idx > -1 && _point_idx < MAX_POINT 
          ;  
}

const double* OpticalRecorder::GetRecord(const double* _pp, int _trk_idx, int _point_idx)  // static
{
    const double* rec = Valid(_trk_idx, _point_idx) ? _pp + 16*MAX_POINT*_trk_idx + 16*_point_idx  : nullptr ; 
    return rec ;     
}




const double* OpticalRecorder::getRecord(int _point_idx) const
{
    return GetRecord( const_cast<const double*>(pp), trk_idx, _point_idx );  
}

void OpticalRecorder::recordPoint( const G4StepPoint* point )
{
    unsigned flag = point_idx == 0 ? unsigned(TORCH) : PointFlag(point) ;  

    if( flag == NAN_ABORT ) return ; 
 
    writePoint( point, flag ); 
     

    if( trk_idx < 5 && 0 ) std::cout 
        << "OpticalRecorder::recordPoint" 
        << " trk_idx " << trk_idx 
        << " point_idx " << point_idx 
        << " flag [" << Flag(flag) << "]"  
        << "\n"
        << descPoint( point_idx )
        << "\n"
        ;  

    point_idx += 1 ; 
}



double OpticalRecorder::GetSpeed(const double* _pp, int _trk_idx, int _point_idx ) // static
{
    const double* curr = GetRecord( _pp, _trk_idx, _point_idx ); 
    const double* prev = GetRecord( _pp, _trk_idx, _point_idx - 1 ); 

    double dt  = DeltaTime( prev, curr ); 
    double dp  = DeltaPos(  prev, curr ); 
    double speed = dt > 0. && dp > 0. ? dp/dt : -1. ; 

    return speed ; 
}


void OpticalRecorder::GetPointSpeed( std::vector<double>& speeds, const double* _pp, uint64_t q_seqhis, int _point_idx )
{
    for(int i=0 ; i < MAX_PHOTON ; i++)
    {
        int _trk_idx = i ; 
        uint64_t seqhis = GetHistory(_pp, _trk_idx);   
        if(seqhis != q_seqhis) continue ;
        double speed = GetSpeed( _pp, _trk_idx, _point_idx );   
        assert( speed >  -1. ); 
        speeds.push_back(speed);  
    }
}


double OpticalRecorder::getSpeed(int _point_idx) const
{ 
    return GetSpeed( pp, trk_idx, _point_idx ); 
}


std::string OpticalRecorder::descPoint(int _point_idx) const
{
    const double* curr = getRecord( _point_idx ); 
    const double* prev = getRecord( _point_idx - 1 ); 

    double dt  = DeltaTime( prev, curr ); 
    double dp  = DeltaPos(  prev, curr ); 
    double speed = dt > 0. && dp > 0. ? dp/dt : -1. ; 
    double speed2 = getSpeed(_point_idx) ; 
    assert( speed == speed2 ); 


    std::stringstream ss ; 
    ss << Desc( curr, 16 ) << "\n" << " dt " << dt << " dp " << dp << " dp/dt " << speed << "\n" ;  
    std::string str = ss.str(); 
    return str ; 
}


uint64_t OpticalRecorder::GetHistory(const double* _pp,  int _trk_idx)
{
    uint64_t seq = 0 ; 
    for(int i=0 ; i < MAX_POINT ; i++)
    {
        int _point_idx = i ; 
        const double* rec = GetRecord(_pp, _trk_idx, _point_idx ); 
        double wl = Wavelength(rec); 
        if(wl == 0.) break ; 
        unsigned fl = PointFlag(rec); 
        seq |= ( ( FFS(fl) & 0xfull ) << i*4 ) ;  
    }
    return seq ; 
}


std::string OpticalRecorder::DescHistory(const double* _pp,  int _trk_idx) 
{
    std::stringstream ss ; 
    for(int i=0 ; i < MAX_POINT ; i++)
    {
        int _point_idx = i ; 
        const double* rec = GetRecord(_pp, _trk_idx, _point_idx ); 
        double wl = Wavelength(rec); 
        if(wl == 0.) break ; 
        unsigned fl = PointFlag(rec); 
        ss << Flag(fl) << " " ; 
    }
    std::string str = ss.str(); 
    return str ; 
}


std::string OpticalRecorder::descHistory(int _trk_idx) const
{
    return DescHistory(pp, _trk_idx) ; 
}
uint64_t OpticalRecorder::getHistory(int _trk_idx) const
{
    uint64_t seq = GetHistory(pp, _trk_idx) ; 
    uint64_t seq1 = getSeq(_trk_idx) ; 
    assert( seq == seq1 ); 
    return seq ; 
}






