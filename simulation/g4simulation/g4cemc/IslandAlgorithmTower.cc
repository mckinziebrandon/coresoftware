
#include "IslandAlgorithmTower.h"
// Local includes.
#include "RawTower.h"
#include "RawTowerGeomContainer.h"
#include "RawTowerContainer.h"

// C++ includes.
#include <iostream>


int IslandAlgorithmTower::maxPhiBin = -10;
int IslandAlgorithmTower::maxEtaBin = -10;

void IslandAlgorithmTower::setCenter(RawTowerGeomContainer* towerGeom) {
    etaCenter = towerGeom->get_etacenter(bineta);
    phiCenter = towerGeom->get_phicenter(binphi);
}

/* ----------------------------------------
   IslandAlgorithmTower(RawTower*) constructor.
   ---------------------------------------- */
IslandAlgorithmTower::IslandAlgorithmTower(RawTower *rt) : id(-1) {
    bineta = rt->get_bineta();
    binphi = rt->get_binphi();
    energy = rt->get_energy();
    id = rt->get_id();
}


// note: true for diagonally adjacent
bool IslandAlgorithmTower::isAdjacent(IslandAlgorithmTower &tower) {
    if (bineta - 1 <= tower.getBinEta() && tower.getBinEta()<=bineta+1) {
        if(binphi-1<=tower.getBinPhi() && tower.getBinPhi()<=binphi+1) {
            return true;
        } else if(((tower.getBinPhi() == maxPhiBin-1) && (binphi == 0)) || 
                ((tower.getBinPhi() == 0) && (binphi == maxPhiBin-1))) {
            return true;
        }
    }
    return false;
}

// Comparison first on bineta if not equal, else on binphi.
bool operator<(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b) {
    if (a.getBinEta() != b.getBinEta()) {
        return a.getBinEta() < b.getBinEta();
    }
    return a.getBinPhi() < b.getBinPhi();
}

bool operator==(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b) {
    return a.getID() == b.getID();
}

RawTower* IslandAlgorithmTower::GetRawTower(IslandAlgorithmTower towerHelper, RawTowerContainer* towers) {
    int iPhi = towerHelper.getBinPhi();
    int iEta = towerHelper.getBinEta();
    // Ensure ids match. TODO: not really sure what this means?
    int eid = (int) RawTowerDefs::encode_towerid(towers->getCalorimeterID(), iEta, iPhi);
    if (towerHelper.getID() != eid) ExitOnIDMismatch(towerHelper.getID(), eid);
    return towers->getTower(iEta, iPhi);
}

void IslandAlgorithmTower::ExitOnIDMismatch(int id1, int id2) {
    std::cout <<__PRETTY_FUNCTION__
        << " - Fatal Error - id mismatch. internal: " 
        << id1
        << ", towercontainer: " 
        << id2
        << std::endl;
    exit(1);
}


