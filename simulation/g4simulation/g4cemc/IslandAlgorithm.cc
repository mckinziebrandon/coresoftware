

#include "IslandAlgorithm.h"

// Local includes.
#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"

// C++ includes.
#include <iostream>

namespace IslandAlgorithm {
    using std::cout;
    using std::endl;


    /* ------------------------------------------------------------------------------------------ *
       1.   The island algorithm starts by a search for seeds. Seeds are defined as                
            crystals with an energy above a certain threshold on transverse energy.                
     * ------------------------------------------------------------------------------------------ */
    std::list<IslandAlgorithmTower> GetSeedTowers(RTContainer* _towers, RTGeomContainer* _towerGeom, float _threshold) {

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
        seedTowers.sort(lessEnergy);
        // Todo: this could _definitely_ be optimized. 
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
        seedTowers.sort(moreEnergy);

        PrintSeeds(seedTowers);
        return seedTowers;
    }

    // Make 5x5 clusters centered on each seed. 
    TowerMap GetSimpleClusters(std::list<IslandAlgorithmTower> seedTowers, 
                                RTContainer*        _towers, 
                                RTGeomContainer*    _towerGeom) {

        TowerMap clusteredTowers;
        int clusterID = 0;
        foreach (IslandAlgorithmTower& seed, seedTowers) {
            clusteredTowers.insert(std::make_pair(clusterID, seed));
            _PrintTowerMsg(seed, clusteredTowers.size(), "SEED");
            int currBinPhi   = seed.getBinPhi();
            int currBinEta   = seed.getBinEta();

            std::vector<int> deltaBins;
            deltaBins.push_back(2);
            deltaBins.push_back(1);
            deltaBins.push_back(0);
            deltaBins.push_back(-1);
            deltaBins.push_back(-2);
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

        int clusterID = 0;
        TowerMap clusteredTowers;
        foreach (IslandAlgorithmTower& seed, seedTowers) {
            // Begin by inserting the seed tower, which defines a cluster.
            clusteredTowers.insert(std::make_pair(clusterID, seed));
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
        RawTower* nextTower = _towers->getTower(currBinEta, _movePhi(direction, currBinPhi));

        // Keep doing this until energy increase or hole.
        if (nextTower && currEnergy > nextTower->get_energy()) {
            IslandAlgorithmTower towerHelper(nextTower);
            towerHelper.setCenter(_towerGeom);
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

        if (_moveEta(direction, currBinEta) != -1) {
            RawTower* nextTower = _towers->getTower(currBinEta, currBinPhi);
            if (nextTower && currEnergy > nextTower->get_energy()) {
                IslandAlgorithmTower towerHelper(nextTower);
                towerHelper.setCenter(_towerGeom);
                clusteredTowers.insert(std::make_pair(clusterID, towerHelper));
                _PrintTowerMsg(towerHelper, clusteredTowers.size(), "ETA");
                _SearchPhi("north", towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchPhi("south", towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
                _SearchEta(direction, towerHelper, clusterID, _towers, clusteredTowers, _towerGeom);
            }
        }
    }

    int _movePhi(std::string direction, int& currBinPhi) { 
        if (direction == "north")       currBinPhi = (currBinPhi != IslandAlgorithmTower::getMaxPhiBin()) ? currBinPhi+1 : 1;
        else if (direction == "south")  currBinPhi = (currBinPhi != 1) ? currBinPhi-1 :
            IslandAlgorithmTower::getMaxPhiBin();
        return currBinPhi;
    }

    int _moveEta(std::string direction, int& currBinEta) { 
        if (direction == "west")        currBinEta = (currBinEta != 1) ? currBinEta-1 : -1;
        else if (direction == "east")   currBinEta = (currBinEta != IslandAlgorithmTower::getMaxEtaBin()) ? currBinEta+1 : -1;
        return currBinEta;
    }

    // A simple comparator that orders IslandAlgorithmTowers in order of INCREASING energy.
    bool lessEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2) { 
        return tower1.getEnergy() < tower2.getEnergy(); 
    }

    // A simple comparator that orders IslandAlgorithmTowers in order of DECREASING energy.
    bool moreEnergy(IslandAlgorithmTower tower1, IslandAlgorithmTower tower2) {
        return !lessEnergy(tower1, tower2);
    }

    // Essentially a 'ToString' method for a list of seed towers.
    void PrintSeeds(std::list<IslandAlgorithmTower>& seeds) {
        foreach (IslandAlgorithmTower& seed, seeds) {
            cout << "seed (energy, eta, phi) = (" 
                 << seed.getEnergy() << ", "
                 << seed.getEtaCenter() << ", "
                 << seed.getPhiCenter() << ")" 
                 << endl;
        }
    }

    void _PrintTowerMsg(IslandAlgorithmTower tower, int index, const char* phiOrEta) {
        cout << index << ". [" << phiOrEta << "] "
             << "Inserted towerID "  << tower.getID()     << endl
             << "\tEnergy="       << tower.getEnergy() << "; "
             << "Eta="           << tower.getEtaCenter() << "; "
             << "Phi="           << tower.getPhiCenter() 
             << endl;
    }

}
