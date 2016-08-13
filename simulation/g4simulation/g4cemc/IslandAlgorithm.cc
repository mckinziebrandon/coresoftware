#include "IslandAlgorithm.h"

// Local includes.
#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"

// C++ includes.
#include <iostream>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace IslandAlgorithm {
    using std::cout;
    using std::endl;

    // Set of towerIDs to ensure towers only associated with one cluster.
    std::set<int> usedTowerIDs;
    int numDoubleCountsAvoided = 0; 

    /* ------------------------------------------------------------------------------------------ *
       1.   The island algorithm starts by a search for seeds. Seeds are defined as                
            crystals with an energy above a certain threshold on transverse energy.                
     * ------------------------------------------------------------------------------------------ */
    std::list<IslandAlgorithmTower> GetSeedTowers(RTContainer* _towers, RTGeomContainer* _towerGeom, float _threshold) {

        typedef std::pair<const unsigned int, RawTower*> RawTowerPair;
        // Collect all towers above threshold.
        std::list<IslandAlgorithmTower> seedTowers;
        foreach (RawTowerPair& rawTowerPair, _towers->getTowers()) {
            if (rawTowerPair.second->get_energy() > _threshold) {
                IslandAlgorithmTower towerHelper(rawTowerPair.second);
                towerHelper.setID(rawTowerPair.first);
                towerHelper.setCenter(_towerGeom);
                seedTowers.push_back(towerHelper);
            }
        }

        // Find towers with higher-energy adjacent towers.(aka bad seeds)
        std::set<IslandAlgorithmTower> badSeeds;
        seedTowers.sort(_LessEnergy);
        // TODO: this could _definitely_ be optimized. 
        foreach (IslandAlgorithmTower& tower1, seedTowers) {
            foreach (IslandAlgorithmTower& tower2, seedTowers) {
                if (tower1.isAdjacent(tower2) && tower1.getEnergy() < tower2.getEnergy()) {
                    badSeeds.insert(tower1);
                }
            }
        }

        // Remove those seeds that are adjacent to higher energy ones.
        foreach (const IslandAlgorithmTower badSeed, badSeeds) {
            seedTowers.remove(badSeed);
        }

        // Order from hi-to-lo energy.
        seedTowers.sort(_MoreEnergy);

        _PrintSeeds(seedTowers);
        return seedTowers;
    }

    // Make 5x5 clusters centered on each seed. 
    TowerMap GetSimpleClusters(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom) {

        TowerMap clusteredTowers;
        int clusterID = 0;
        foreach (IslandAlgorithmTower& seed, seedTowers) {
            _PrintTowerMsg(seed, clusteredTowers.size(), "SEED");
            clusteredTowers.insert(std::make_pair(clusterID, seed));
            int currBinPhi   = seed.getBinPhi();
            int currBinEta   = seed.getBinEta();

            int deltaVals[] = {2, 1, 0, -1, -2};
            std::vector<int> deltaBins (deltaVals, deltaVals + sizeof(deltaVals) / sizeof(int));
            foreach (int& deltaPhiBin, deltaBins) {
                foreach (int& deltaEtaBin, deltaBins) {
                    int binPhi = currBinPhi + deltaPhiBin;
                    int binEta = currBinEta + deltaEtaBin;
                    if (deltaPhiBin == 0 && deltaEtaBin == 0) continue;
                    RawTower* rawTower = _towers->getTower(binEta, binPhi);
                    if (rawTower) {
                        IslandAlgorithmTower towerHelper(rawTower);
                        towerHelper.setCenter(_towerGeom);
                        clusteredTowers.insert(std::make_pair(clusterID, towerHelper));
                    }
                }
            }
            clusterID++;
        }
        return clusteredTowers;
    }

    /* ------------------------------------------------------------------------------------------ *
       2.   Starting from the most energetic seed, the algorithm collects crystals belonging to 
            a certain cluster. Moves both directions in phi, collecting all towers until it sees
            a rise in energy, or a hole. The algorithm then steps in eta and makes another phi
            search. The eta-steps are stopped when a rise in energy, or a hole, is encountered.
            When in direction in eta is completed, the algorithm goes back to the seed position 
            and works in the other eta direction. All towers are makred as belonging to that one 
            cluster and can't be subsequently used to seed another. 
     * ------------------------------------------------------------------------------------------ */
    TowerMap GetClusteredTowers(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom) {

        // Initialize/declare required new objects.
        int clusterID = 0;
        numDoubleCountsAvoided = 0; 
        usedTowerIDs.clear();
        TowerMap clusteredTowers;

        foreach (IslandAlgorithmTower& seed, seedTowers) {
            // Begin by inserting the seed tower, which defines a cluster.
            clusteredTowers.insert(std::make_pair(clusterID, seed));
            usedTowerIDs.insert(seed.getID());

            _PrintTowerMsg(seed, clusteredTowers.size(), "SEED");
            _SearchPhi("north", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchPhi("south", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchEta("west", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            _SearchEta("east", seed, clusterID, _towers, clusteredTowers, _towerGeom);
            clusterID++;
        }

        return clusteredTowers;
    }


    /* -------------------------------------------------------------------------
       _SearchPhi
       ------------------------------------------------------------------------- */
    void _SearchPhi(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom) {

        // Get current tower info.
        float currEnergy = currentTower.getEnergy();
        int currBinPhi   = currentTower.getBinPhi();
        int currBinEta   = currentTower.getBinEta();

        // Get next tower info to decide if we should add it.
        RawTower* nextTower = _towers->getTower(currBinEta, _MovePhi(direction, currBinPhi, _towerGeom));

        // Keep doing this until energy increase or hole.
        if (nextTower && currEnergy > nextTower->get_energy()) {
            IslandAlgorithmTower towerHelper(nextTower);
            towerHelper.setCenter(_towerGeom);

            // Terminate search if encounter tower that has already been clustered.
            if (_TowerAlreadyClustered(towerHelper)) return;
            // Otherwise, mark it as clustered and insert it into the clusteredTowers map.
            usedTowerIDs.insert(towerHelper.getID());
            clusteredTowers.insert(std::make_pair(clusterID, towerHelper));

            _PrintTowerMsg(towerHelper, clusteredTowers.size(), "PHI");
            _SearchPhi(direction, towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
        }
    }


    /* -------------------------------------------------------------------------
       _SearchEta
       ------------------------------------------------------------------------- */
    void _SearchEta(std::string direction,       IslandAlgorithmTower           currentTower,
                    int&        clusterID,       RTContainer*       _towers,    
                    TowerMap&   clusteredTowers, RTGeomContainer*   _towerGeom) {

        // Get next tower info to decide if we should add it.
        int currBinEta   = currentTower.getBinEta();
        int currBinPhi   = currentTower.getBinPhi();
        float currEnergy = currentTower.getEnergy();

        if (_MoveEta(direction, currBinEta, _towerGeom) != -1) {
            RawTower* nextTower = _towers->getTower(currBinEta, currBinPhi);
            if (nextTower && currEnergy > nextTower->get_energy()) {
                IslandAlgorithmTower towerHelper(nextTower);
                towerHelper.setCenter(_towerGeom);

                // Terminate search if encounter tower that has already been clustered.
                if (_TowerAlreadyClustered(towerHelper)) return;
                // Otherwise, mark it as clustered and insert it into the clusteredTowers map.
                usedTowerIDs.insert(towerHelper.getID());
                clusteredTowers.insert(std::make_pair(clusterID, towerHelper));

                _PrintTowerMsg(towerHelper, clusteredTowers.size(), "ETA");
                _SearchPhi("north", towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchPhi("south", towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchEta(direction, towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
            }
        }
    }

    int _MovePhi(std::string direction, int& currBinPhi, RTGeomContainer* _towerGeom) { 
        if (direction == "north")       currBinPhi = (currBinPhi != _towerGeom->get_phibins()) ? currBinPhi+1 : 1;
        else if (direction == "south")  currBinPhi = (currBinPhi != 1) ? currBinPhi-1 :
            _towerGeom->get_phibins();
        return currBinPhi;
    }

    int _MoveEta(std::string direction, int& currBinEta, RTGeomContainer* _towerGeom) { 
        if (direction == "west")        currBinEta = (currBinEta != 1) ? currBinEta-1 : -1;
        else if (direction == "east")   currBinEta = (currBinEta != _towerGeom->get_etabins()) ? currBinEta+1 : -1;
        return currBinEta;
    }


    // A simple comparator that orders IslandAlgorithmTowers in order of INCREASING energy.
    bool _LessEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2) { 
        return tower1.getEnergy() < tower2.getEnergy(); 
    }

    // A simple comparator that orders IslandAlgorithmTowers in order of DECREASING energy.
    bool _MoreEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2) {
        return !_LessEnergy(tower1, tower2);
    }

    // Essentially a 'ToString' method for a list of seed towers.
    void _PrintSeeds(std::list<IslandAlgorithmTower>& seeds) {
        foreach (IslandAlgorithmTower& seed, seeds) {
            cout << "seed (energy, eta, phi) = (" 
                 << seed.getEnergy() << ", "
                 << seed.getEtaCenter() << ", "
                 << seed.getPhiCenter() << ")" 
                 << endl;
        }
    }


    void _PrintTowerMsg(IslandAlgorithmTower tower, int index, const char* towerType) {
        cout << index << ". [" << towerType << "] "
             << "Inserted towerID "  << tower.getID()     << endl
             << "\tEnergy="       << tower.getEnergy() << "; "
             << "Eta="           << tower.getEtaCenter() << "; "
             << "Phi="           << tower.getPhiCenter() 
             << endl;
    }

    bool _TowerAlreadyClustered(IslandAlgorithmTower towerHelper) {
        bool res = (usedTowerIDs.find(towerHelper.getID()) != usedTowerIDs.end());
        if (res) _PrintDebugMsg(towerHelper);
        return res;
    }

    void _PrintDebugMsg(IslandAlgorithmTower towerHelper) {  
        numDoubleCountsAvoided++;
        cout << "\n -------------------------------------------------------- "  << endl
             << " DOUBLE-COUNT "    << numDoubleCountsAvoided << " AVOIDED. "   << endl
             << "Tower ID:\t"       << towerHelper.getID()              << endl
             << "Tower Eta:\t"      << towerHelper.getEtaCenter()       << endl
             << "Tower Phi:\t"      << towerHelper.getPhiCenter()       << endl;
    }

}
