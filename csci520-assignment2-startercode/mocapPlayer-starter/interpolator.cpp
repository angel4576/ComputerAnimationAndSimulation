#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "types.h"

#include <iostream>
#include <vector>

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for (int frame = 1; frame <= N; frame++)
    {
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
    }

    startKeyframe = endKeyframe;

  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));


}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
    double roll = angles[0] * M_PI / 180.0; // X
    double pitch = angles[1] * M_PI / 180.0; // Y
    double yaw = angles[2] * M_PI / 180.0;   // Z

    double cx = cos(roll), sx = sin(roll);  
    double cy = cos(pitch), sy = sin(pitch); 
    double cz = cos(yaw), sz = sin(yaw);

    // R = Rz * Ry * Rx (row-major) 
    R[0] = cz * cy;    R[1] = cz * sy * sx - sz * cx;    R[2] = cz * sy * cx + sz * sx;
    R[3] = sz * cy;    R[4] = sz * sy * sx + cz * cx;    R[5] = sz * sy * cx - cz * sx;
    R[6] = -sy;        R[7] = cy * sx;                   R[8] = cy * cx;

}

vector Bezier3D(vector p0, vector p1, vector p2, vector p3, double t)
{
    double u = 1 - t;
    vector result(u * u * u * p0.x() + 3 * u * u * t * p1.x() + 3 * u * t * t * p2.x() + t * t * t * p3.x()
        , u * u * u * p0.y() + 3 * u * u * t * p1.y() + 3 * u * t * t * p2.y() + t * t * t * p3.y()
        , u * u * u * p0.z() + 3 * u * u * t * p1.z() + 3 * u * t * t * p2.z() + t * t * t * p3.z());
    
    return result;
}

vector NormalizeAngle(vector angle) {
    for (int i = 0; i < 3; i++) {
        while (angle.p[i] > 180) angle.p[i] -= 360;
        while (angle.p[i] < -180) angle.p[i] += 360;
    }
    return angle;
}

std::vector<vector> ComputerBezierControlPointsEuler(vector prev, vector start, vector end, vector next)
{
    std::vector<vector> controlPts(2);

    // calculate control points
    vector eExtrapolated = (start - prev) * 2.0 + start;
    vector aBar = (eExtrapolated + end) * 0.5;
    vector a = start + (aBar - start) / 3.0;

    vector eExtrapolated1 = (end - start) * 2.0 + end;
    vector a1Bar = (eExtrapolated1 + next) * 0.5;
    vector a1 = end + (a1Bar - end) / 3.0;
    vector b1 = end + (a1Bar - end) / -3.0;

    vector controlPoints[] = { a, b1 };
    controlPts[0] = a;
    controlPts[1] = b1;

    return controlPts;
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
    // total frame number
    int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

    int startKeyframe = 0;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;

        Posture* startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture* endPosture = pInputMotion->GetPosture(endKeyframe);

        int prevKeyframe = startKeyframe - (N + 1);
        int nextKeyframe = endKeyframe + N + 1;
        if (prevKeyframe < 0)
        {
            prevKeyframe = 0;
        }

        if (nextKeyframe >= inputLength)
        {
            nextKeyframe = inputLength - 1;
        }

        Posture* prevPosture = pInputMotion->GetPosture(prevKeyframe);
        Posture* nextPosture = pInputMotion->GetPosture(nextKeyframe);

        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);


        // interpolate in between
        for (int frame = 1; frame <= N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N + 1);

            // interpolate root position
            
            // calculate control points
            std::vector<vector> controlPts = ComputerBezierControlPointsEuler(prevPosture->root_pos, startPosture->root_pos, endPosture->root_pos, nextPosture->root_pos);     
           
            interpolatedPosture.root_pos = DeCasteljauEuler(t, startPosture->root_pos, controlPts[0], controlPts[1], endPosture->root_pos);
            
            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
            {               
                // calculate control points
                std::vector<vector> controlPts = ComputerBezierControlPointsEuler(prevPosture->bone_rotation[bone], startPosture->bone_rotation[bone], endPosture->bone_rotation[bone], nextPosture->bone_rotation[bone]);
         
                interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, startPosture->bone_rotation[bone], controlPts[0], controlPts[1], endPosture->bone_rotation[bone]);
                
            }

            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }

        startKeyframe = endKeyframe;

    }

    for (int frame = startKeyframe + 1; frame < inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
    int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

    int startKeyframe = 0;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;

        Posture* startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture* endPosture = pInputMotion->GetPosture(endKeyframe);

        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);

        // interpolate in between
        for (int frame = 1; frame <= N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N + 1);

            // interpolate root position
            interpolatedPosture.root_pos = startPosture->root_pos * (1 - t) + endPosture->root_pos * t;

            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
            {
                Quaternion<double> qStart, qEnd;
                Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
                Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd);

                // test
                /*std::cout << "---- before ----" << std::endl;
                std::cout << startPosture->bone_rotation[bone].p[0] << ", "
                    << startPosture->bone_rotation[bone].p[1] << ", "
                    << startPosture->bone_rotation[bone].p[2] << std::endl;

                Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart);
                double test[3];
                Quaternion2Euler(qStart, test);

                std::cout << "---- after ----" << std::endl;
                std::cout << test[0] << ", "
                    << test[1] << ", "
                    << test[2] << std::endl;

                return;*/

                Quaternion<double> qInterpolated;
                // Slerp
                qInterpolated = Slerp(t, qStart, qEnd);
                  
                Quaternion2Euler(qInterpolated, interpolatedPosture.bone_rotation[bone].p);

            }

            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }

        startKeyframe = endKeyframe;

    }

    for (int frame = startKeyframe + 1; frame < inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

}

Quaternion<double> Bisect(Quaternion<double> p, Quaternion<double> q)
{
    Quaternion<double> sum = p + q;
    double norm = sum.Norm();

    Quaternion<double> result = sum / norm;
    return result;
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
    // total frame number
    int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

    int startKeyframe = 0;
    while (startKeyframe + N + 1 < inputLength)
    {
        int endKeyframe = startKeyframe + N + 1;

        Posture* startPosture = pInputMotion->GetPosture(startKeyframe);
        Posture* endPosture = pInputMotion->GetPosture(endKeyframe);

        int prevKeyframe = startKeyframe - (N + 1);
        int nextKeyframe = endKeyframe + N + 1;
        if (prevKeyframe < 0)
        {
            prevKeyframe = 0;
        }

        if (nextKeyframe >= inputLength)
        {
            nextKeyframe = inputLength - 1;
        }

        Posture* prevPosture = pInputMotion->GetPosture(prevKeyframe);
        Posture* nextPosture = pInputMotion->GetPosture(nextKeyframe);

        // copy start and end keyframe
        pOutputMotion->SetPosture(startKeyframe, *startPosture);
        pOutputMotion->SetPosture(endKeyframe, *endPosture);

        // interpolate in between
        for (int frame = 1; frame <= N; frame++)
        {
            Posture interpolatedPosture;
            double t = 1.0 * frame / (N + 1);

            // interpolate root position
            std::vector<vector> controlPts = ComputerBezierControlPointsEuler(prevPosture->root_pos, startPosture->root_pos, endPosture->root_pos, nextPosture->root_pos);
            
            // Bezier Euler
            interpolatedPosture.root_pos = DeCasteljauEuler(t, startPosture->root_pos, controlPts[0], controlPts[1], endPosture->root_pos);

            // interpolate bone rotations
            for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
            {
                Quaternion<double> qPrev, qStart, qEnd, qNext;
                Euler2Quaternion(prevPosture->bone_rotation[bone].p, qPrev); // qn-1
                Euler2Quaternion(startPosture->bone_rotation[bone].p, qStart); // qn
                Euler2Quaternion(endPosture->bone_rotation[bone].p, qEnd); // qn+1
                Euler2Quaternion(nextPosture->bone_rotation[bone].p, qNext); // qn+2

                // calculate control points
                // Quaternion<double> qExtrapolated = Slerp(2.0, qPrev, qStart);
                Quaternion<double> qExtrapolated = Double(qPrev, qStart);
                Quaternion<double> aBar = Slerp(0.5, qExtrapolated, qEnd);  
                Quaternion<double> a = Slerp(1.0 / 3, qStart, aBar); 
                // Quaternion<double> b = Slerp(-1.0 / 3, qStart, aBar); 

                Quaternion<double> qExtrapolated1 = Double(qStart, qEnd);
                Quaternion<double> a1Bar = Slerp(0.5, qExtrapolated1, qNext);
                Quaternion<double> a1 = Slerp(1.0 / 3, qEnd, a1Bar);
                Quaternion<double> b1 = Slerp(-1.0 / 3, qEnd, a1Bar);


                // Shoemake paper calculation
                //Quaternion<double> dPrevStart = Double(qPrev, qStart);
                //Quaternion<double> an = Bisect(dPrevStart, qEnd); // a_n = Bisect(Double(q_n-1, qn), q_n+1)
                //
                //Quaternion<double> dStartEnd = Double(qStart, qEnd); 
                //Quaternion<double> an1 = Bisect(dStartEnd, qNext);  
                //Quaternion<double> bn1 = Double(an1, qEnd); // b_n+1  

                //an = Slerp(1.0 / 3, qStart, an);  
                //bn1 = Slerp(1.0 / 3, qEnd, bn1);


                Quaternion<double> qInterpolated = DeCasteljauQuaternion(t, qStart, a, b1, qEnd);
             
                Quaternion2Euler(qInterpolated, interpolatedPosture.bone_rotation[bone].p);

            }

            pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
        }

        startKeyframe = endKeyframe;

    }

    for (int frame = startKeyframe + 1; frame < inputLength; frame++)
        pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
    double R[9];
    Euler2Rotation(angles, R);
    
    q = Quaternion<double>::Matrix2Quaternion(R);
    q.Normalize();
    
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
    double R[9];
    q.Quaternion2Matrix(R);

    Rotation2Euler(R, angles);
     
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_)
{
    // qStart dot qEnd
    double startDotEnd = qStart.Gets() * qEnd_.Gets() + qStart.Getx() * qEnd_.Getx()
                        + qStart.Gety() * qEnd_.Gety() + qStart.Getz() * qEnd_.Getz();

    if (startDotEnd < 0.0) // theta > 90
    {
        qEnd_ = -1 * qEnd_; // make sure go to shortest path
        startDotEnd = -startDotEnd;
    }

    // clamp
    if (startDotEnd < -1.0)
    {
        startDotEnd = -1.0;
    }
    else if(startDotEnd > 1.0)
    {
        startDotEnd = 1.0;
    }

    double theta = acos(startDotEnd); 
    double sinTheta = sin(theta);

    if (sinTheta < 1e-5) {
        return qStart;
    }

    // nan check
    /*if (std::isnan(theta))
    {
        std::cerr << "ERROR: theta result is NaN!" << std::endl; 
    }

    if (std::isnan(sinTheta))
    {
        std::cerr << "ERROR: sinTheta result is NaN!" << std::endl; 
    }*/

    double w1 = sin((1 - t) * theta) / sinTheta; 
    double w2 = sin(t * theta) / sinTheta;

    // nan check
    /*if (std::isnan(w1))
    {
        std::cerr << "ERROR: w1 result is NaN!" << std::endl;
    }

    if (std::isnan(w2))
    {
        std::cerr << "ERROR: w1 result is NaN!" << std::endl;
    }*/

    Quaternion<double> result;
    result = w1 * qStart + w2 * qEnd_;

    return result;
}



Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
    // Double(p, q) = 2(p¡¤q)q ¨C p
    double pDotq = p.Gets() * q.Gets() + p.Getx() * q.Getx()
        + p.Gety() * q.Gety() + p.Getz() * q.Getz();

    Quaternion<double> result;
    result = 2 * pDotq * q - p;

    return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
    vector q0 = p0 * (1 - t) + p1 * t;
    vector q1 = p1 * (1 - t) + p2 * t;
    vector q2 = p2 * (1 - t) + p3 * t;

    vector r0 = q0 * (1 - t) + q1 * t;
    vector r1 = q1 * (1 - t) + q2 * t;

    vector result = r0 * (1 - t) + r1 * t;

    return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
    Quaternion<double> q0 = Slerp(t, p0, p1); 
    Quaternion<double> q1 = Slerp(t, p1, p2); 
    Quaternion<double> q2 = Slerp(t, p2, p3);  
    
    Quaternion<double> r0 = Slerp(t, q0, q1); 
    Quaternion<double> r1 = Slerp(t, q1, q2); 

    Quaternion<double> result = Slerp(t, r0, r1);
    return result;
}

