#include "IK.h"
#include "FK.h"
#include "minivectorTemplate.h"
#include <Eigen/Dense>
#include <adolc/adolc.h>
#include <cassert>
#if defined(_WIN32) || defined(WIN32)
  #ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
  #endif
#endif
#include <math.h>
using namespace std;

// CSCI 520 Computer Animation and Simulation
// Jernej Barbic and Yijing Li

namespace
{

// Converts degrees to radians.
template<typename real>
inline real deg2rad(real deg) { return deg * M_PI / 180.0; }

template<typename real>
Mat3<real> Euler2Rotation(const real angle[3], RotateOrder order)
{
  Mat3<real> RX = Mat3<real>::getElementRotationMatrix(0, deg2rad(angle[0]));
  Mat3<real> RY = Mat3<real>::getElementRotationMatrix(1, deg2rad(angle[1]));
  Mat3<real> RZ = Mat3<real>::getElementRotationMatrix(2, deg2rad(angle[2]));

  switch(order)
  {
    case RotateOrder::XYZ:
      return RZ * RY * RX;
    case RotateOrder::YZX:
      return RX * RZ * RY;
    case RotateOrder::ZXY:
      return RY * RX * RZ;
    case RotateOrder::XZY:
      return RY * RZ * RX;
    case RotateOrder::YXZ:
      return RZ * RX * RY;
    case RotateOrder::ZYX:
      return RX * RY * RZ;
  }
  assert(0);
}

// Performs forward kinematics, using the provided "fk" class.
// This is the function whose Jacobian matrix will be computed using adolc.
// numIKJoints and IKJointIDs specify which joints serve as handles for IK:
//   IKJointIDs is an array of integers of length "numIKJoints"
// Input: numIKJoints, IKJointIDs, fk, eulerAngles (of all joints)
// Output: handlePositions (world-coordinate positions of all the IK joints; length is 3 * numIKJoints)
template<typename real>
void forwardKinematicsFunction(
    int numIKJoints, const int * IKJointIDs, const FK & fk,
    const std::vector<real> & eulerAngles, std::vector<real> & handlePositions)
{
  // Students should implement this.
  // The implementation of this function is very similar to function computeLocalAndGlobalTransforms in the FK class.
  // The recommended approach is to first implement FK::computeLocalAndGlobalTransforms.
  // Then, implement the same algorithm into this function. To do so,
  // you can use fk.getJointUpdateOrder(), fk.getJointRestTranslation(), and fk.getJointRotateOrder() functions.
  // Also useful is the multiplyAffineTransform4ds function in minivectorTemplate.h .
  // It would be in principle possible to unify this "forwardKinematicsFunction" and FK::computeLocalAndGlobalTransforms(),
  // so that code is only written once. We considered this; but it is actually not easily doable.
  // If you find a good approach, feel free to document it in the README file, for extra credit.

    int numJoints = fk.getNumJoints();
    std::vector<Mat3<real>> globalRotations(numJoints);
    std::vector<Vec3<real>> globalTranslations(numJoints);
    
    for (int i = 0; i < numJoints; i++)
    {
        int jointId = fk.getJointUpdateOrder(i);
        
        // Local transform
        real jointEulerAngle[3] = {
            eulerAngles[3 * jointId + 0], // eulerAngle.size = 3 * numOfJoints
            eulerAngles[3 * jointId + 1],
            eulerAngles[3 * jointId + 2]
        };
        Mat3<real> rotation = Euler2Rotation(jointEulerAngle, fk.getJointRotateOrder(jointId));
        
        // Joint orientation
        const Vec3d& jointOrient = fk.getJointOrient(jointId);
        real orientEulerAngle[3] = {
            static_cast<real>(jointOrient[0]), // eulerAngle.size = 3 * numOfJoints
            static_cast<real>(jointOrient[1]),
            static_cast<real>(jointOrient[2]) 
        };
        Mat3<real> orientRotation = Euler2Rotation(orientEulerAngle, fk.getJointRotateOrder(jointId));

        Mat3<real> localRotation = orientRotation * rotation;
        
        const Vec3d& jointTrans = fk.getJointRestTranslation(jointId); 
        Vec3<real> localTranslation(
            static_cast<real>(jointTrans[0]),
            static_cast<real>(jointTrans[1]),
            static_cast<real>(jointTrans[2])
        );
        
        // multiplyAffineTransform4ds
        // RigidTransform4d localTransform =  RigidTransform4d(orientRotation * rotation, fk.getJointRestTranslation(jointId));


        // Global transform
        int parentId = fk.getJointParent(jointId);
        if (parentId >= 0) // not root
        {
            // globalTransform = parentGlobalTransform * localTransform
            // globalRotations[parentId] 
            multiplyAffineTransform4ds(globalRotations[parentId], globalTranslations[parentId]
            , localRotation, localTranslation
            , globalRotations[jointId], globalTranslations[jointId]); // Rout, tout
        }
        else
        {
            globalRotations[jointId] = localRotation;
            globalTranslations[jointId] = localTranslation;
        }

    }

    // Get world positions of IK joints from global translations
    // handlePositions.resize(3 * numIKJoints);
    for (int i = 0; i < numIKJoints; i++)
    {
        int IKJointId = IKJointIDs[i];
        const Vec3<real>& position = globalTranslations[IKJointId];

        handlePositions[3 * i + 0] = position[0];
        handlePositions[3 * i + 1] = position[1]; 
        handlePositions[3 * i + 2] = position[2]; 
    }


}

} // end anonymous namespaces

IK::IK(int numIKJoints, const int * IKJointIDs, FK * inputFK, int adolc_tagID)
{
  this->numIKJoints = numIKJoints;
  this->IKJointIDs = IKJointIDs;
  this->fk = inputFK;
  this->adolc_tagID = adolc_tagID;

  FKInputDim = fk->getNumJoints() * 3;
  FKOutputDim = numIKJoints * 3;

  train_adolc();
}

void IK::train_adolc()
{
  // Students should implement this.
  // Here, you should setup adol_c:
  //   Define adol_c inputs and outputs. 
  //   Use the "forwardKinematicsFunction" as the function that will be computed by adol_c.
  //   This will later make it possible for you to compute the gradient of this function in IK::doIK
  //   (in other words, compute the "Jacobian matrix" J).
  // See ADOLCExample.cpp .

    trace_on(adolc_tagID);

    std::vector<adouble> x(FKInputDim); // define input of function f
    std::vector<adouble> y(FKOutputDim);

    for (int i = 0; i < FKInputDim; i++)
    {
        x[i] <<= 0.0;
    }

    // The computation of f (FK)
    forwardKinematicsFunction(numIKJoints, IKJointIDs, *fk, x, y); 
    
    std::vector<double> output(FKOutputDim);
    for (int i = 0; i < FKOutputDim; i++)
    {
        y[i] >>= output[i];
    }

    trace_off();

}

void IK::doIK(const Vec3d * targetHandlePositions, Vec3d * jointEulerAngles)
{
  // You may find the following helpful:
  int numJoints = fk->getNumJoints(); // Note that is NOT the same as numIKJoints!

  // Students should implement this.
  // Use adolc to evalute the forwardKinematicsFunction and its gradient (Jacobian). It was trained in train_adolc().
  // Specifically, use ::function, and ::jacobian .
  // See ADOLCExample.cpp .
  //
  // Use it implement the Tikhonov IK method (or the pseudoinverse method for extra credit).
  // Note that at entry, "jointEulerAngles" contains the input Euler angles. 
  // Upon exit, jointEulerAngles should contain the new Euler angles.

  int n = FKInputDim; // input dimension 
  int m = FKOutputDim; // output dimension

  std::vector<double> theta(FKInputDim);
  for (int i = 0; i < numJoints; i++)
  {
      theta[3 * i + 0] = jointEulerAngles[i][0];
      theta[3 * i + 1] = jointEulerAngles[i][1];
      theta[3 * i + 2] = jointEulerAngles[i][2];
  }

  // Compute handle positions of IK joints by FK (f(¦È))
  std::vector<double> currentHandlePosition(FKOutputDim);
  ::function(adolc_tagID, FKOutputDim, FKInputDim, theta.data(), currentHandlePosition.data());

  // Compute ¦¤b (diff between target pos and current pos/ b - f(¦È))
  std::vector<double> posDiff(FKOutputDim);
  for (int i = 0; i < numIKJoints; i++)
  {
      posDiff[3 * i + 0] =  targetHandlePositions[i][0] - currentHandlePosition[3 * i + 0];
      posDiff[3 * i + 1] =  targetHandlePositions[i][1] - currentHandlePosition[3 * i + 1];
      posDiff[3 * i + 2] =  targetHandlePositions[i][2] - currentHandlePosition[3 * i + 2];
  }

  // Compute J 
  std::vector<double> jacobianMatrix(FKOutputDim * FKInputDim); // m x n
  std::vector<double*> jacobianMatrixEachRow(FKOutputDim); // pointer array where each pointer points to one row of the jacobian matrix
  for (int i = 0; i < FKOutputDim; i++)
  {
      jacobianMatrixEachRow[i] = &jacobianMatrix[i * FKInputDim];
  }
  ::jacobian(adolc_tagID, FKOutputDim, FKInputDim, theta.data(), jacobianMatrixEachRow.data()); 

  // Solve linear system 
  Eigen::MatrixXd J(FKOutputDim, FKInputDim); // m x n J matrix
  for (int i = 0; i < m; i++)
  {
      for (int j = 0; j < n; j++)
      {
          J(i, j) = jacobianMatrix[i * n + j];
      }

  }

  // (J^T * J + a * I) * ¦¤¦È = J^T * ¦¤b
  // Left-hand side
  double alpha = 1e-3;
  Eigen::MatrixXd A = J.transpose()* J + alpha * Eigen::MatrixXd::Identity(n, n); 

  // Right-hand side
  Eigen::VectorXd delta_b(m, 1); 
  for (int i = 0; i < m; i++) 
  {
      delta_b[i] = posDiff[i];
  }
    
  Eigen::VectorXd rhs = J.transpose() * delta_b; 

  // Solve for x in A x = b
  Eigen::VectorXd delta_theta = A.ldlt().solve(rhs);
  
  for (int i = 0; i < n; i++) 
  {
      theta[i] += delta_theta[i];
      jointEulerAngles[i / 3][i % 3] = theta[i];
  }

}

