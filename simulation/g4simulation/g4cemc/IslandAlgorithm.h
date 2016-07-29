#ifndef ISLANDALGORITHM_H_
#define ISLANDALGORITHM_H_

#include "IslandAlgorithmTower.h"
#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"

#include <map>
#include <set>



class RawTower;
class RawTowerGeomContainer;
class RawTower;

namespace IslandAlgorithm {
    using std::cout;
    using std::endl;

    // Jin: note if you want to make these class swappable for other uses, please use template function
    typedef RawTowerContainer               RTContainer;
    typedef RawTowerGeomContainer           RTGeomContainer;

    typedef std::multimap<int, IslandAlgorithmTower>    TowerMap;
    typedef std::pair<const int, IslandAlgorithmTower>  TowerPair;
    typedef std::pair<const unsigned int, RawTower*> RawTowerPair;

    // The two main workhorse functions of the island algorithm. 
    std::list<IslandAlgorithmTower> GetSeedTowers(RTContainer* _towers, 
                                      RTGeomContainer* _towerGeom, 
                                      float _threshold=0.);

    // Make 5x5 clusters centered on each seed. 
    TowerMap GetSimpleClusters(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom);

    TowerMap GetClusteredTowers(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer* _towers, 
                                RTGeomContainer* _towerGeom);

    // Collect towers in specified phi/eta direction for specified cluster.
    void _SearchPhi(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);
    void _SearchEta(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom);

    // Advance phi/eta depending on current location.
    int _movePhi(std::string direction, int& currBinPhi);
    int _moveEta(std::string direction, int& currBinEta);
    // Comparators for sorting.
    bool lessEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2);
    bool moreEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2);
    // Basic print function for energy, etacenter, phicenter.
    void PrintSeeds(std::list<IslandAlgorithmTower>& seeds);
    void _PrintTowerMsg(IslandAlgorithmTower tower, int index, const char* phiOrEta);


}

#endif
