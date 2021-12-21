//*********************************************************
//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
//*********************************************************

//**********************************************************************************************
//
// VolumetricPrimitives.hlsli
//
// Ray marching of Metaballs (aka "Blobs").
// More info here: https://www.scratchapixel.com/lessons/advanced-rendering/rendering-distance-fields/blobbies
//
//**********************************************************************************************

#ifndef VOLUMETRICPRIMITIVESLIBRARY_H
#define VOLUMETRICPRIMITIVESLIBRARY_H


#include "RaytracingShaderHelper.hlsli"

struct Metaball
{
    float3 center;
    float  radius;
};

//**test**zsym
struct Group
{
    int head;
    int tail;
    float near;
    float far;
};


// Calculate a magnitude of an influence from a Metaball charge.
// Return metaball potential range: <0,1>
// mbRadius - largest possible area of metaball contribution - AKA its bounding sphere.
float CalculateMetaballPotential(in float3 position, in Metaball blob, out float distance)
{
    distance = length(position - blob.center);
    
    if (distance <= blob.radius)
    {
        float d = distance;

        // Quintic polynomial field function.
        // The advantage of this polynomial is having smooth second derivative. Not having a smooth
        // second derivative may result in a sharp and visually unpleasant normal vector jump.
        // The field function should return 1 at distance 0 from a center, and 1 at radius distance,
        // but this one gives f(0) = 0, f(radius) = 1, so we use the distance to radius instead.
        d = blob.radius - d;

        float r = blob.radius;
        return 6 * (d*d*d*d*d) / (r*r*r*r*r)
            - 15 * (d*d*d*d) / (r*r*r*r)
            + 10 * (d*d*d) / (r*r*r);
    }
    return 0;
}

// Calculate field potential from all active metaballs.
float CalculateMetaballsPotential(in float3 position, in Metaball blobs[N_METABALLS], in UINT nActiveMetaballs)
{
    float sumFieldPotential = 0;
#if USE_DYNAMIC_LOOPS 
    for (UINT j = 0; j < nActiveMetaballs; j++)
#else
    for (UINT j = 0; j < N_METABALLS; j++)
#endif
    {
        float dummy;
        sumFieldPotential += CalculateMetaballPotential(position, blobs[j], dummy);
    }
    return sumFieldPotential;
}

// Calculate a normal via central differences.
float3 CalculateMetaballsNormal(in float3 position, in Metaball blobs[N_METABALLS], in UINT nActiveMetaballs)
{
    float e = 0.5773 * 0.00001;
    return normalize(float3(
        CalculateMetaballsPotential(position + float3(-e, 0, 0), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(e, 0, 0), blobs, nActiveMetaballs),
        CalculateMetaballsPotential(position + float3(0, -e, 0), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(0, e, 0), blobs, nActiveMetaballs),
        CalculateMetaballsPotential(position + float3(0, 0, -e), blobs, nActiveMetaballs) -
        CalculateMetaballsPotential(position + float3(0, 0, e), blobs, nActiveMetaballs)));
}

void InitializeAnimatedMetaballs(out Metaball blobs[N_METABALLS], in float elapsedTime, in float cycleDuration)
{


    // Metaball centers at t0 and t1 key frames.
#if N_METABALLS == 5
    float3 keyFrameCenters[N_METABALLS][2] =
    {
        { float3(-0.7, 0, 0),float3(0.7,0, 0) },
        { float3(0.7 , 0, 0), float3(-0.7, 0, 0) },
        { float3(0, -0.7, 0),float3(0, 0.7, 0) },
        { float3(0, 0.7, 0), float3(0, -0.7, 0) },
        { float3(0, 0, 0),   float3(0, 0, 0) }
    };
    // Metaball field radii of max influence
    float radii[N_METABALLS] = { 0.35, 0.35, 0.35, 0.35, 0.25 };

    //#if N_METABALLS == 10
    //float3 keyFrameCenters[N_METABALLS][2] =
    //{
    //    { float3(-0.5, 0, 0),float3(0, 0.5, 0) },
    //    { float3(0.7, 0, 0), float3(-0.7, 0, 0) },
    //    { float3(0, -0.7, 0),float3(0, 0.3, 0) },
    //    { float3(0, 0.1, 0), float3(0, -0.6, 0) },
    //    { float3(0, 0, 0),   float3(0, 0, 1) },
    //    { float3(0, 0, 0.3),float3(0.6, 0.2, 0) },
    //    { float3(0, 0, 0.3), float3(0.2, 0, 0) },
    //    { float3(0, 0.3, 0),float3(-0.5, 0.3, 0) },
    //    { float3(0, -0.1, 0.2), float3(0, 0.6, 0) },
    //    { float3(1, 0, 0),   float3(0, 0, 0.1) }
    //};
    //// Metaball field radii of max influence
    //float radii[N_METABALLS] = { 0.35, 0.35, 0.35, 0.35, 0.25, 0.35, 0.35, 0.35, 0.35, 0.25 };

#else
    float3 keyFrameCenters[N_METABALLS][2] =
    {
        //{ float3(-0.3, -0.3, -0.4),float3(0.3,-0.3,-0.0) },
        //{ float3(0.0, -0.2, 0.5), float3(0.0, 0.4, 0.5) },
        //{ float3(0.4,0.4, 0.4), float3(-0.4, 0.2, -0.4) }
        { float3(0.0, 0.0, 0.0), float3(1.0, 0.0, 0.0) },
        { float3(0.0, 0.0, 0.0), float3(0.0, 1.0, 0.0) },
        { float3(0.0, 0.0, 0.0), float3(0.0, 0.0, 1.0) }
    };
    // Metaball field radii of max influence
    float radii[N_METABALLS] = { 0.40, 0.50, 0.60 };
#endif

    // Calculate animated metaball center positions.
    float  tAnimate = CalculateAnimationInterpolant(elapsedTime, cycleDuration);
    for (UINT j = 0; j < N_METABALLS; j++)
    {
        blobs[j].center = lerp(keyFrameCenters[j][0], keyFrameCenters[j][1], tAnimate);
        blobs[j].radius = radii[j];
    }
}

// Find all metaballs that ray intersects.
// The passed in array is sorted to the first nActiveMetaballs.
//**rayt**void FindIntersectingMetaballs(in Ray ray, out float tmin, out float tmax, inout Metaball blobs[N_METABALLS], out UINT nActiveMetaballs, \
//    out float blobsTmin[N_METABALLS], out float blobsTmax[N_METABALLS])
void FindIntersectingMetaballs(in Ray ray, out float tmin, out float tmax, inout Metaball blobs[N_METABALLS], \
    inout float centerDistance[N_METABALLS], out UINT nActiveMetaballs)
{
    // Find the entry and exit points for all metaball bounding spheres combined.
    tmin = INFINITY;
    tmax = -INFINITY;

    nActiveMetaballs = 0;
    for (UINT i = 0; i < N_METABALLS; i++)
    {
        float _thit, _tmax;
        if (RaySolidSphereIntersectionTest(ray, _thit, _tmax, blobs[i].center, blobs[i].radius))
        {
            tmin = min(_thit, tmin);
            tmax = max(_tmax, tmax);
#if LIMIT_TO_ACTIVE_METABALLS
            blobs[nActiveMetaballs] = blobs[i];
            centerDistance[nActiveMetaballs] = centerDistance[i];
            //blobsTmin[nActiveMetaballs] = _thit;
            //blobsTmax[nActiveMetaballs] = _tmax;
            nActiveMetaballs++;
#else
            nActiveMetaballs = N_METABALLS;
            //blobsTmin[i] = _thit;
            //blobsTmax[i] = _tmax;
#endif
        }
    }
    tmin = max(tmin, RayTMin());
    tmax = min(tmax, RayTCurrent());
}

//**test** by zsym
//**old**
//void sortMetaballsByDistance0(in Ray ray, in UINT nActiveMetaballs, inout Metaball blobs[N_METABALLS], \
//    inout float centerDistance[N_METABALLS], inout float blobsTmin[N_METABALLS], inout float blobsTmax[N_METABALLS])
//{
//#if USE_DYNAMIC_LOOPS
//    for (UINT k = 0; k < nActiveMetaballs; k++)
//#else
//    for (UINT k = 0; k < N_METABALLS; k++)
//#endif
//    {
//        float3 originToMetaballCenter = blobs[k].center - ray.origin;
//        centerDistance[k] = dot(originToMetaballCenter, normalize(ray.direction));
//    }
//
//    // Bubble sort
//#if USE_DYNAMIC_LOOPS
//    for (UINT i = 0; i < nActiveMetaballs - 1; i++)
//    {
//        for (UINT j = 0; j < nActiveMetaballs - 1 - i; j++)
//#else
//    for (UINT i = 0; i < N_METABALLS - 1; i++)
//    {
//        for (UINT j = 0; j < N_METABALLS - 1 - i; j++)
//#endif
//        {
//            if (centerDistance[j] > centerDistance[j + 1])
//            {
//                /*float tmpD = centerDistance[j];
//                centerDistance[j] = centerDistance[j + 1];
//                centerDistance[j + 1] = tmpD;*/
//
//                swap(centerDistance[j], centerDistance[j + 1]);
//
//                Metaball tmpB = blobs[j];
//                blobs[j] = blobs[j + 1];
//                blobs[j + 1] = tmpB;
//
//                float tmpTmin = blobsTmin[j];
//                blobsTmin[j] = blobsTmin[j + 1];
//                blobsTmin[j + 1] = tmpTmin;
//
//                float tmpTmax = blobsTmax[j];
//                blobsTmax[j] = blobsTmax[j + 1];
//                blobsTmax[j + 1] = tmpTmax;
//            }
//        }
//    }
//
//
//    // Bubble sort, single loop version, slower than the dual loop version
////#if USE_DYNAMIC_LOOPS
////    UINT bound = nActiveMetaballs - 1;
////#else
////    UINT bound = N_METABALLS - 1;
////#endif
////
////    for (UINT j = 0; j < bound; j++)
////    {
////        if (centerDistance[j] > centerDistance[j + 1])
////        {
////            /*float tmpD = centerDistance[j];
////            centerDistance[j] = centerDistance[j + 1];
////            centerDistance[j + 1] = tmpD;*/
////
////            swap(centerDistance[j], centerDistance[j + 1]);
////
////            Metaball tmpB = blobs[j];
////            blobs[j] = blobs[j + 1];
////            blobs[j + 1] = tmpB;
////
////            float tmpTmin = blobsTmin[j];
////            blobsTmin[j] = blobsTmin[j + 1];
////            blobsTmin[j + 1] = tmpTmin;
////
////            float tmpTmax = blobsTmax[j];
////            blobsTmax[j] = blobsTmax[j + 1];
////            blobsTmax[j + 1] = tmpTmax;
////        }
////
////        if (j == (bound - 1))
////        {
////            --j;
////            --bound;
////        }
////    }
//}



//**test** by zsym, this function uses STATIC LOOPS, is called before FindIntersectingMetaballs().
void SortMetaballsByDistance(in Ray ray, inout Metaball blobs[N_METABALLS], out float centerDistance[N_METABALLS])
{
    for (UINT k = 0; k < N_METABALLS; k++)
    {
        float3 originToMetaballCenter = blobs[k].center - ray.origin;
        centerDistance[k] = dot(originToMetaballCenter, normalize(ray.direction));
    }

    // Bubble sort, dual-loop version
    for (UINT i = 0; i < N_METABALLS - 1; i++)
    {
        for (UINT j = 0; j < N_METABALLS - 1 - i; j++)
        {
            if (centerDistance[j] > centerDistance[j + 1])
            {
                swap(centerDistance[j], centerDistance[j + 1]);

                Metaball tmpB = blobs[j];
                blobs[j] = blobs[j + 1];
                blobs[j + 1] = tmpB;
            }
        }
    }
}

//void divideGroups(in centerDistance[N_METABALLS], in float blobsTmin[N_METABALLS], in float blobsTmax[N_METABALLS, \
//    out Group groups[N_GROUPS]) {}

void DivideGroups(in float centerDistance[N_METABALLS], in Metaball blobs[N_METABALLS], \
    out Group groups[N_GROUPS], in UINT nActiveMetaballs, out UINT groupNum)
{
    // `groupNum` always equals the index of the group to add new blob to
    groupNum = -1;

#if USE_DYNAMIC_LOOPS
    for (UINT i = 0; i < nActiveMetaballs; i++)
#else
    for (UINT i = 0; i < N_METABALLS; i++)
#endif
    {
        // Whether the new blob should be added to current group
        if ((groupNum >= 0) && (((centerDistance[i] - blobs[i].radius) < groups[groupNum].far) || (groupNum == (N_GROUPS - 1))))
        {
            // Append the new blob into current group
            groups[groupNum].tail = i;
            groups[groupNum].far = max(groups[groupNum].far, (centerDistance[i] + blobs[i].radius));

            if ((centerDistance[i] - blobs[i].radius) < groups[groupNum].near)
            { 
                // The new blob refreshes the group's `near`
                groups[groupNum].near = centerDistance[i] - blobs[i].radius;

                while (groupNum > 0)
                { 
                    // Whether the new `near` leads to overlapping of former groups
                    if (groups[groupNum - 1].far > groups[groupNum].near)
                    { 
                        // Merge two groups
                        groups[groupNum - 1].tail = groups[groupNum].tail;
                        groups[groupNum - 1].far = groups[groupNum].far;
                        // Decrement group counter
                        --groupNum;
                    }
                }
            }
        }

        // Create a new group
        else
        {
            ++groupNum;
            groups[groupNum].head = i;
            groups[groupNum].tail = i;
            groups[groupNum].near = centerDistance[i] - blobs[i].radius;
            groups[groupNum].far = centerDistance[i] + blobs[i].radius;
        }
            
    }
}

// Test if a ray with RayFlags and segment <RayTMin(), RayTCurrent()> intersects metaball field.
// The test sphere traces through the metaball field until it hits a threshold isosurface. 
bool RayMetaballsIntersectionTest(in Ray ray, out float thit, out ProceduralPrimitiveAttributes attr, in float elapsedTime)
{
    Metaball blobs[N_METABALLS];
    InitializeAnimatedMetaballs(blobs, elapsedTime, 12.0f);

    float centerDistance[N_METABALLS];
    SortMetaballsByDistance(ray, blobs, centerDistance);
    
    float tmin, tmax;   // Ray extents to first and last metaball intersections.
    UINT nActiveMetaballs = 0;  // Number of metaballs's that the ray intersects.
    
    FindIntersectingMetaballs(ray, tmin, tmax, blobs, centerDistance, nActiveMetaballs);

    Group groups[N_GROUPS];
    UINT groupNum = 0;
    DivideGroups(centerDistance, blobs, groups, nActiveMetaballs, groupNum);
    //**old**
    //float centerDistance[N_METABALLS];
    //sortMetaballsByDistance0(ray, nActiveMetaballs, blobs, centerDistance, blobsTmin, blobsTmax);

    UINT MAX_STEPS = 128;
    float t = tmin;
    float minTStep = (tmax - tmin) / (MAX_STEPS / 1);
    UINT iStep = 0;

    while (iStep++ < MAX_STEPS)
    {
        float3 position = ray.origin + t * ray.direction;
        float fieldPotentials[N_METABALLS];    // Field potentials for each metaball.
        float sumFieldPotential = 0;           // Sum of all metaball field potentials.
            
        // Calculate field potentials from all metaballs.
#if USE_DYNAMIC_LOOPS
        for (UINT j = 0; j < nActiveMetaballs; j++)
#else
        for (UINT j = 0; j < N_METABALLS; j++)
#endif
        {
            float distance;
            fieldPotentials[j] = CalculateMetaballPotential(position, blobs[j], distance);
            sumFieldPotential += fieldPotentials[j];
         }

        // Field potential threshold defining the isosurface.
        // Threshold - valid range is (0, 1>, the larger the threshold the smaller the blob.
        const float Threshold = 0.25f;

        // Have we crossed the isosurface?
        if (sumFieldPotential >= Threshold)
        {
            float3 normal = CalculateMetaballsNormal(position, blobs, nActiveMetaballs);
            if (IsAValidHit(ray, t, normal))
            {
                thit = t;
                attr.normal = normal;
                return true;
            }
        }
        t += minTStep;
    }

    return false;
}

#endif // VOLUMETRICPRIMITIVESLIBRARY_H