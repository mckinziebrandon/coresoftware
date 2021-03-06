#ifndef PHG4CEmcTestBeamSubsystem_h
#define PHG4CEmcTestBeamSubsystem_h

#include "g4main/PHG4Subsystem.h"

#include <Geant4/G4Types.hh>

class PHCompositeNode;
class PHG4CEmcTestBeamDetector;
class PHG4CEmcTestBeamSteppingAction;
class PHG4EventAction;

class PHG4CEmcTestBeamSubsystem: public PHG4Subsystem
{

  public:

  //! constructor
  PHG4CEmcTestBeamSubsystem( const std::string &name = "BLOCK", const int layer = 0 );

  //! destructor
  virtual ~PHG4CEmcTestBeamSubsystem( void )
  {}

  //! init
  /*!
  creates the detector_ object and place it on the node tree, under "DETECTORS" node (or whatever)
  reates the stepping action and place it on the node tree, under "ACTIONS" node
  creates relevant hit nodes that will be populated by the stepping action and stored in the output DST
  */
  int Init(PHCompositeNode *);

  //! event processing
  /*!
  get all relevant nodes from top nodes (namely hit list)
  and pass that to the stepping action
  */
  int process_event(PHCompositeNode *);

  //! accessors (reimplemented)
  virtual PHG4Detector* GetDetector( void ) const;
  virtual PHG4SteppingAction* GetSteppingAction( void ) const;

  void SetSize(const G4double sizex, const G4double sizey, const G4double sizez)
     {dimension[0] = sizex; dimension[1] = sizey; dimension[2] = sizez;}
  void SetPlaceZ(const G4double dbl);
  void SetPlace(const G4double place_x, const G4double place_y, const G4double place_z);
  void SetXRot(const G4double dbl);
  void SetYRot(const G4double dbl);
  void SetZRot(const G4double dbl);
  PHG4EventAction* GetEventAction() const {return eventAction_;}
  void SetActive(const int i = 1) {active = i;}
  void SetAbsorberActive(const int i = 1) {absorberactive = i;}
  void SuperDetector(const std::string &name) {superdetector = name;}
  const std::string SuperDetector() {return superdetector;}

  void BlackHole(const int i=1) {blackhole = i;}

  private:

  //! detector geometry
  /*! defives from PHG4Detector */
  PHG4CEmcTestBeamDetector* detector_;

  //! particle tracking "stepping" action
  /*! derives from PHG4SteppingActions */
  PHG4CEmcTestBeamSteppingAction* steppingAction_;
  PHG4EventAction *eventAction_;
  G4double dimension[3];
  G4double place_in_x;
  G4double place_in_y;
  G4double place_in_z;
  G4double rot_in_x;
  G4double rot_in_y;
  G4double rot_in_z;

  int active;
  int absorberactive;
  int layer;
  int blackhole;
  std::string detector_type;
  std::string superdetector;
};

#endif
