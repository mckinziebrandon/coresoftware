#ifndef ISLANDALGORITHM_H_
#define ISLANDALGORITHM_H_

#include "IslandAlgorithmTower.h"

#include <map>
#include <set>
#include <list>


class RawTower;
class RawTowerGeomContainer;
class RawTowerContainer;


namespace IslandAlgorithm {
    using std::cout;
    using std::endl;

    // Jin: note if you want to make these class swappable for other uses, please use template function
    typedef RawTowerContainer               RTContainer;
    typedef RawTowerGeomContainer           RTGeomContainer;

    typedef std::multimap<int, IslandAlgorithmTower>    TowerMap;
    typedef std::pair<const int, IslandAlgorithmTower>  TowerPair;
    //typedef std::pair<const unsigned int, RawTower*> RawTowerPair;


    // The two main workhorse functions of the island algorithm. 
    std::list<IslandAlgorithmTower> GetSeedTowers(RTContainer* _towers, 
                                      RTGeomContainer* _towerGeom, 
                                      float _threshold=0.);


    TowerMap GetClusteredTowers(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer* _towers, 
                                RTGeomContainer* _towerGeom);

    // Alternative to GetClusteredTowers; Make 5x5 clusters centered on each seed. 
    TowerMap GetSimpleClusters(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom);


    // -------------------------------------------------------------------------
    // Helper Functions
    // -------------------------------------------------------------------------

    // Collect towers in specified phi/eta direction for specified cluster.
    void _SearchPhi(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);
    void _SearchEta(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);

    // Advance phi/eta depending on current location.
    //int _MovePhi(std::string direction, int& currBinPhi, );
    int _MovePhi(std::string direction, int& currBinPhi, RTGeomContainer* _towerGeom);
    int _MoveEta(std::string direction, int& currBinEta, RTGeomContainer* _towerGeom);
    //int _MoveEta(std::string direction, int& currBinEta);

    // Comparators for sorting.
    bool _LessEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2);
    bool _MoreEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2);

    // Tell whether or not a tower has been assigned to a cluster this event.
    bool _TowerAlreadyClustered(IslandAlgorithmTower towerHelper);

    // Basic print function for energy, etacenter, phicenter.
    void _PrintSeeds(std::list<IslandAlgorithmTower>& seeds);
    void _PrintTowerMsg(IslandAlgorithmTower tower, int index, const char* towerType);
    void _PrintDebugMsg(IslandAlgorithmTower towerHelper);

}

#endif
