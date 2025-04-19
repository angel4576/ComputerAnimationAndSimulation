/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include "input.h"

#include <iostream>


/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    // initialize acceleration array
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                a[i][j][k].x = 0.0;
                a[i][j][k].y = 0.0;
                a[i][j][k].z = 0.0;
            }
    
    double structRestLength = 1.0 / 7.0;
    double bendRestLength = 2.0 * structRestLength;

    for (int i = 0; i <= 7; i++)
        for (int j = 0; j <= 7; j++)
            for (int k = 0; k <= 7; k++)
            {                
                // Calculate structural springs 
                if (i < 7) 
                {
                    calculateSpring(jello, a, i, j, k, i + 1, j, k, structRestLength);
                    

                    // Calculate bend spring 
                    if (i < 6)
                    {
                        calculateSpring(jello, a, i, j, k, i + 2, j, k, bendRestLength);
                        

                    }
                    
                }
               
                if (j < 7)
                {
                    calculateSpring(jello, a, i, j, k, i, j + 1, k, structRestLength);
                    

                    // Calculate bend spring 
                    if (j < 6)
                    {
                        calculateSpring(jello, a, i, j, k, i, j + 2, k, bendRestLength);                       

                    }

                }
                
                if (k < 7)
                {
                    calculateSpring(jello, a, i, j, k, i, j, k + 1, structRestLength);
                    

                    // Calculate bend spring 
                    if (k < 6)
                    {
                        calculateSpring(jello, a, i, j, k, i, j, k + 2, bendRestLength);
                        

                    }
                }
                
                // Calculate shear springs
                calculateShearSpring(jello, a, i, j, k);                

                // Calculate force field
                point externalForce = interpolateForce(jello, i, j, k);
                a[i][j][k].x += externalForce.x / jello->mass;
                a[i][j][k].y += externalForce.y / jello->mass;
                a[i][j][k].z += externalForce.z / jello->mass;
                                
                /*
                    Equal force controlled by mouse
               
                */ 
                if (g_isDragging)
                {
                    point uniformDragForce = {dragForce[0], dragForce[1], dragForce[2]};
                    double length = sqrt(uniformDragForce.x * uniformDragForce.x + uniformDragForce.y * uniformDragForce.y + uniformDragForce.z * uniformDragForce.z);
                    if (length > 0)
                    {
                        pNORMALIZE(uniformDragForce);
                    }

                    // clamp
                    uniformDragForce.x = (uniformDragForce.x >= 2.0) ? 2.0 : uniformDragForce.x;
                    uniformDragForce.y = (uniformDragForce.y >= 2.0) ? 2.0 : uniformDragForce.y;
                    uniformDragForce.z = (uniformDragForce.z >= 2.0) ? 2.0 : uniformDragForce.z;
                                       
                    double factor = 1;
                    a[i][j][k].x += uniformDragForce.x * factor / jello->mass; 
                    a[i][j][k].y += uniformDragForce.y * factor / jello->mass;
                    a[i][j][k].z += uniformDragForce.z * factor / jello->mass;

                    /*char buffer[256]; 
                    sprintf(buffer, "Drag: x = %f, y = %f, z = %f\n", uniformDragForce.x, uniformDragForce.y, uniformDragForce.z);
                    OutputDebugStringA(buffer);*/
                }
                
                
                /*
                    Collision
                */

                // aabb collision
                if (isCollideWithBox(jello->p[i][j][k], 4.0))
                {
                    // OutputDebugStringA("Collide\n");
                    point penaltyForce = computePenaltyForce(jello, jello->p[i][j][k], jello->v[i][j][k]);

                    a[i][j][k].x += penaltyForce.x / jello->mass;
                    a[i][j][k].y += penaltyForce.y / jello->mass;
                    a[i][j][k].z += penaltyForce.z / jello->mass;
                }

                // inclined plane collision
                if (isCollideWithPlane(jello->p[i][j][k], jello->a, jello->b, jello->c, jello->d))
                {
                    // OutputDebugStringA("Collide with plane!!\n");
                    point planePenaltyForce = computePlanePenaltyForce(jello, jello->p[i][j][k], jello->v[i][j][k]);
                    a[i][j][k].x += planePenaltyForce.x / jello->mass;
                    a[i][j][k].y += planePenaltyForce.y / jello->mass;
                    a[i][j][k].z += planePenaltyForce.z / jello->mass;
                }

                    
            }

}

// Trilinear Interpolation
point interpolateForce(struct world* jello, int i, int j, int k)
{
    // n x n x n force field (n: resolution)    
    int n = jello->resolution;

    double gridSize = 4.0 / (n - 1);

    double gridX = (jello->p[i][j][k].x + 2.0) / gridSize; // p - (-2(bound)) / grid size 
    double gridY = (jello->p[i][j][k].y + 2.0) / gridSize; 
    double gridZ = (jello->p[i][j][k].z + 2.0) / gridSize; 

    // bottom-left corner of the grid cube
    int i0 = (int)gridX;
    int j0 = (int)gridY; 
    int k0 = (int)gridZ;
    
    // top-right corner
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;

    // avoid array out of bound 
    if (i1 >= n) i1 = i0;
    if (j1 >= n) j1 = j0;
    if (k1 >= n) k1 = k0;

    // get 8 points of the cube   
    point F000 = jello->forceField[i0 * n * n + j0 * n + k0];
    point F100 = jello->forceField[i1 * n * n + j0 * n + k0];
    point F010 = jello->forceField[i0 * n * n + j1 * n + k0];
    point F110 = jello->forceField[i1 * n * n + j1 * n + k0];
    
    point F001 = jello->forceField[i0 * n * n + j0 * n + k1];
    point F101 = jello->forceField[i1 * n * n + j0 * n + k1];
    point F011 = jello->forceField[i0 * n * n + j1 * n + k1];
    point F111 = jello->forceField[i1 * n * n + j1 * n + k1];

    // weight
    double tx = gridX - i0;
    double ty = gridY - j0;
    double tz = gridZ - k0;

    // trilinear interpolation
    // X
    point F00 = { (1 - tx) * F000.x + tx * F100.x, 
        (1 - tx) * F000.y + tx * F100.y, 
        (1 - tx) * F000.z + tx * F100.z };

    point F10 = { (1 - tx) * F010.x + tx * F110.x,
        (1 - tx) * F010.y + tx * F110.y,
        (1 - tx) * F010.z + tx * F110.z };

    point F01 = { (1 - tx) * F001.x + tx * F101.x,
        (1 - tx) * F001.y + tx * F101.y,
        (1 - tx) * F001.z + tx * F101.z };

    point F11 = { (1 - tx) * F011.x + tx * F111.x,
        (1 - tx) * F011.y + tx * F111.y,
        (1 - tx) * F011.z + tx * F111.z };

    // Y
    point F0 = { (1 - ty) * F00.x + ty * F10.x,
        (1 - ty) * F00.y + ty * F10.y,
        (1 - ty) * F00.z + ty * F10.z };

    point F1 = { (1 - ty) * F01.x + ty * F11.x,
        (1 - ty) * F01.y + ty * F11.y,
        (1 - ty) * F01.z + ty * F11.z };

    // Z
    point finalForce = { (1 - tz) * F0.x + tz * F1.x,
        (1 - tz) * F0.y + tz * F1.y,
        (1 - tz) * F0.z + tz * F1.z };

    return finalForce;
}

void calculateSpring(struct world* jello, struct point a[8][8][8]
    , int i, int j, int k, int nextI, int nextJ, int nextK, double restLength)
{
    point totalForce = { 0, 0, 0 };

    point springForce = computeSpringForce(jello->p[i][j][k], jello->p[nextI][nextJ][nextK], jello->kElastic, restLength);
    point dampingForce = computeDampingForce(jello->p[i][j][k], jello->p[nextI][nextJ][nextK]
        , jello->v[i][j][k], jello->v[nextI][nextJ][nextK], jello->dElastic);


    totalForce.x = springForce.x + dampingForce.x;
    totalForce.y = springForce.y + dampingForce.y;
    totalForce.z = springForce.z + dampingForce.z;

    a[i][j][k].x += totalForce.x / jello->mass;
    a[i][j][k].y += totalForce.y / jello->mass;
    a[i][j][k].z += totalForce.z / jello->mass;

    // counter force
    // update neighbor acceleration (opposite direction) 
    a[nextI][nextJ][nextK].x -= totalForce.x / jello->mass;
    a[nextI][nextJ][nextK].y -= totalForce.y / jello->mass;
    a[nextI][nextJ][nextK].z -= totalForce.z / jello->mass;
}

void calculateShearSpring(struct world* jello, struct point a[8][8][8], int i, int j, int k)
{
    int shearOffset[4][3] = {
        {1, 1, 0},
        {0, 1, 1},
        {1, 0, 1},
        {1, 1, 1}
    };

    int forwardShearOffset[12][3] = {
        {1, 1, 0}, {-1, 1, 0},  // xy

        {1, 0, 1}, {-1, 0, 1},  // xz
        {1, 0, -1}, {-1, 0, -1},

        {0, 1, 1}, {0, 1, -1},  // yz

        {1, 1, 1}, {-1, 1, 1},
        {1, 1, -1}, {-1, 1, -1}

    };

    int completeShearOffset[16][3] = {
        {1, 1, 0}, {-1, 1, 0}, {1, -1, 0}, {-1, -1, 0}, // xy face diagonal
        {1, 0, 1}, {-1, 0, 1}, {1, 0, -1}, {-1, 0, -1}, // xz 
        {0, 1, 1}, {0, -1, 1}, {0, 1, -1}, {0, -1, -1}, // yz
        {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}   // body diagonal 

    };

    for (int s = 0; s < 12; s++)
    {
        int di = forwardShearOffset[s][0];
        int dj = forwardShearOffset[s][1];
        int dk = forwardShearOffset[s][2];

        // bound check
        if (i + di >= 0 && i + di < 8
            && j + dj >= 0 && j + dj < 8
            && k + dk >= 0 && k + dk < 8)
        {
            // calculate diagnal rest length
            double restLength = sqrt(di * di + dj * dj + dk * dk) / 7.0;

            point shearSpringForce = computeSpringForce(jello->p[i][j][k], jello->p[i + di][j + dj][k + dk]
                , jello->kElastic, restLength);

            point shearDampingForce = computeDampingForce(jello->p[i][j][k], jello->p[i + di][j + dj][k + dk]
                , jello->v[i][j][k], jello->v[i + di][j + dj][k + dk]
                , jello->dElastic);

            point shearTotalForce = { 0, 0, 0 };
            shearTotalForce.x = shearSpringForce.x + shearDampingForce.x;
            shearTotalForce.y = shearSpringForce.y + shearDampingForce.y;
            shearTotalForce.z = shearSpringForce.z + shearDampingForce.z;

            /*char buffer[256];
            sprintf(buffer, "Shear: x = %f, y = %f, z = %f\n", shearTotalForce.x, shearTotalForce.y, shearTotalForce.z);
            OutputDebugStringA(buffer);*/

            a[i][j][k].x += shearTotalForce.x / jello->mass;
            a[i][j][k].y += shearTotalForce.y / jello->mass;
            a[i][j][k].z += shearTotalForce.z / jello->mass;

            // update diagonal neighbor
            a[i + di][j + dj][k + dk].x -= shearTotalForce.x / jello->mass;
            a[i + di][j + dj][k + dk].y -= shearTotalForce.y / jello->mass;
            a[i + di][j + dj][k + dk].z -= shearTotalForce.z / jello->mass;
        }


    }
}

point computeSpringForce(point pA, point pB, double kHook, double restLength)
{
    // vector from p2 to p1
    point L;
    L.x = pA.x - pB.x;
    L.y = pA.y - pB.y;
    L.z = pA.z - pB.z;
    
    double vecLength = sqrt(L.x * L.x + L.y * L.y + L.z * L.z); 
    if (vecLength == 0)
    {
        //vecLength = 1;
    }

    // calculate unit vector L / |L|
    double length;
    pNORMALIZE(L);

    // double restLength = 1.0 / 7.0;
    double mag = -kHook * (vecLength - restLength);
    
    point springForce;
    springForce.x = mag * L.x;
    springForce.y = mag * L.y; 
    springForce.z = mag * L.z;

    return springForce;
}

point computeDampingForce(point pA, point pB, point vA, point vB, double kDamp) 
{
    // vector from pB to pA
    point L;
    L.x = pA.x - pB.x; 
    L.y = pA.y - pB.y; 
    L.z = pA.z - pB.z; 

    double lengthSqr = L.x * L.x + L.y * L.y + L.z * L.z;
    if (lengthSqr == 0)
    {
        //lengthSqr = 1;
    }
    
    point relativeV; // vA - vB
    relativeV.x = vA.x - vB.x;
    relativeV.y = vA.y - vB.y;
    relativeV.z = vA.z - vB.z;

    double vDotL = relativeV.x * L.x + relativeV.y * L.y + relativeV.z * L.z;
    double factor = -kDamp * vDotL / lengthSqr;

    point dampingForce;
    dampingForce.x = factor * L.x;
    dampingForce.y = factor * L.y;
    dampingForce.z = factor * L.z;


    return dampingForce;
}

point computePenaltyForce(struct world* jello, point p, point velocity) // collision response
{
    double penetrationX = 0.0 , penetrationY = 0.0, penetrationZ = 0.0;
    if (p.x > 2) // collide with +x
    {
        penetrationX = p.x - 2;
    }
    else if (p.x < -2)
    {
        penetrationX = -2 - p.x;
    }
   
    // Y
    if (p.y > 2)
    {
        penetrationY = p.y - 2;
    }
    else if (p.y < -2)
    {
        penetrationY = -2 - p.y;
    }

    // Z
    if (p.z > 2)
    {
        penetrationZ = p.z - 2;
    }
    else if (p.z < -2)
    {
        penetrationZ = -2 - p.z;
    }

    /* 
        decide which face is collided
        if collide with multiple faces, choose min penetration one
    */
    double minPenetration = penetrationX;
    int axis = 0; // x: 0, y: 1, z: 2
    if (penetrationY > 0 && (penetrationY < minPenetration || minPenetration == 0))
    {
        minPenetration = penetrationY;
        axis = 1;
    }

    if (penetrationZ > 0 && (penetrationZ < minPenetration || minPenetration == 0))
    {
        minPenetration = penetrationZ;
        axis = 2;
    }

    point contactP = p;
    if (axis == 0) // x
    {
        contactP.x = (p.x > 2) ? 2 : (p.x < -2) ? -2 : p.x; 
    }
    else if (axis == 1)
    {
        contactP.y = (p.y > 2) ? 2 : (p.y < -2) ? -2 : p.y;
    }
    else if (axis == 2)
    {
        contactP.z = (p.z > 2) ? 2 : (p.z < -2) ? -2 : p.z;    
    }


    point totalF = { 0, 0, 0 };
    point springF = computeSpringForce(p, contactP, jello->kCollision, 0.0); // pA: p  
    point dampingF = computeDampingForce(p, contactP, velocity, { 0, 0, 0 }, jello->dCollision);

    totalF.x = springF.x + dampingF.x;
    totalF.y = springF.y + dampingF.y;
    totalF.z = springF.z + dampingF.z;
        
    return totalF;
}

point computePlanePenaltyForce(struct world* jello, point p, point velocity)
{
    double a = jello->a;
    double b = jello->b;
    double c = jello->c;
    double d = jello->d;

    // distance from p to plane
    double nLength = sqrt(a * a + b * b + c * c); 
    double distance = (a * p.x + b * p.y + c * p.z + d) / nLength; // signed

    point contactP; // = { 0, 0, 0 };
    // Pc = P - d * n
    point normal = { a, b, c };
    double length;
    pNORMALIZE(normal);

    contactP.x = p.x - distance * normal.x;
    contactP.y = p.y - distance * normal.y;
    contactP.z = p.z - distance * normal.z;

    // calculate penalty
    point totalF = { 0, 0, 0 };
    point springF = computeSpringForce(p, contactP, jello->kCollision, 0.0); // pA: p   
    point dampingF = computeDampingForce(p, contactP, velocity, { 0, 0, 0 }, jello->dCollision); 

    totalF.x = springF.x + dampingF.x; 
    totalF.y = springF.y + dampingF.y; 
    totalF.z = springF.z + dampingF.z; 

    return totalF;
}

bool isCollideWithBox(point p, double boxWidth)
{
    double halfW = boxWidth / 2.0;
    if (p.x >= -halfW && p.x <= halfW &&
        p.y >= -halfW && p.y <= halfW &&
        p.z >= -halfW && p.z <= halfW)
    {
        return false;
    }
    else
    {
        return true;
    }

}

bool isCollideWithPlane(point p, double a, double b, double c, double d)
{
   double Fp = a * p.x + b * p.y + c * p.z + d;
    
    return Fp < 0; // collide / go to the other side of the plane
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
