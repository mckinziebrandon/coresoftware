// Raw Tower Helper Class.
#ifndef __ISLANDALGORITHMTOWER_H__
#define __ISLANDALGORITHMTOWER_H__

// C++ includes.
#include <iostream>
#include <cmath>
// Local includes.
#include "RawTowerDefs.h"

class RawTower;
class RawTowerGeomContainer;
class RawTowerContainer;


class IslandAlgorithmTower {
    public:
        IslandAlgorithmTower(RawTower *);
        virtual ~IslandAlgorithmTower() {}
        bool isAdjacent(IslandAlgorithmTower &);
        int getID() const                       { return id; }
        int getBinEta() const                   { return bineta; }
        int getBinPhi() const                   { return binphi; }
        float getET() const                     { return energy / std::cosh(etaCenter);}
        float getEnergy() const                 { return energy; }
        float getEtaCenter() const              { return etaCenter; }
        float getPhiCenter() const              { return phiCenter; }
        void setID(const int i)                 { id = i; }
        void setEnergy(float e)                 { energy = e; }
        void setEtaCenter(float eta)            { etaCenter = eta; }
        void setPhiCenter(float phi)            { phiCenter = phi; }
        void setCenter(RawTowerGeomContainer*);
        static int getMaxPhiBin()               { return maxPhiBin; }
        static int getMaxEtaBin()               { return maxEtaBin; }
        static void setMaxPhiBin(const int i)   { maxPhiBin = i; }
        static void setMaxEtaBin(const int i)   { maxEtaBin = i; }
        static RawTower* GetRawTower(IslandAlgorithmTower, RawTowerContainer*);
    protected:
        RawTowerDefs::keytype id;
        int bineta; 
        int binphi;
        float energy;
        float etaCenter; 
        float phiCenter;
        static int maxPhiBin;
        static int maxEtaBin;
        static void ExitOnIDMismatch(int id1, int id2);
};

bool operator<(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b);
bool operator==(const IslandAlgorithmTower& a, const IslandAlgorithmTower& b);

#endif
