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
    UINT head;
    UINT tail;
    float near;
    float far;
    float tmin;
    float tmax;
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

// Calculate field potential from all active metaballs in current group
float CalculateGroupMetaballsPotential(in float3 position, in Metaball blobs[N_METABALLS], UINT head, UINT tail)
{
    float sumFieldPotential = 0;
    for (UINT i = head; i <= tail; ++i)
    {
        float dummy;
        sumFieldPotential += CalculateMetaballPotential(position, blobs[i], dummy);
    }
    return sumFieldPotential;
}

// Calculate a normal via central differences. Only blobs in the group is needed.
float3 CalculateGroupMetaballsNormal(in float3 position, in Metaball blobs[N_METABALLS], UINT head, UINT tail)
{
    float e = 0.5773 * 0.00001;
    return normalize(float3(
        CalculateGroupMetaballsPotential(position + float3(-e, 0, 0), blobs, head, tail) -
        CalculateGroupMetaballsPotential(position + float3(e, 0, 0), blobs, head, tail),
        CalculateGroupMetaballsPotential(position + float3(0, -e, 0), blobs, head, tail) -
        CalculateGroupMetaballsPotential(position + float3(0, e, 0), blobs, head, tail),
        CalculateGroupMetaballsPotential(position + float3(0, 0, -e), blobs, head, tail) -
        CalculateGroupMetaballsPotential(position + float3(0, 0, e), blobs, head, tail)));
}

void InitializeAnimatedMetaballs(out Metaball blobs[N_METABALLS], in float elapsedTime, in float cycleDuration)
{


    // Metaball centers at t0 and t1 key frames.
//#if N_METABALLS == 5
//    float3 keyFrameCenters[N_METABALLS][2] =
//    {
//        { float3(-0.7, 0, 0),float3(0.7,0, 0) },
//        { float3(0.7 , 0, 0), float3(-0.7, 0, 0) },
//        { float3(0, -0.7, 0),float3(0, 0.7, 0) },
//        { float3(0, 0.7, 0), float3(0, -0.7, 0) },
//        { float3(0, 0, 0),   float3(0, 0, 0) }
//    };
//    // Metaball field radii of max influence
//    float radii[N_METABALLS] = { 0.35, 0.35, 0.35, 0.35, 0.25 };
//
#if N_METABALLS == 5
    float3 keyFrameCenters[N_METABALLS][2] =
    {
        { float3(-0.7, 0, 0), float3(0.7,0, 0) },
        { float3(0.7 , 0, 0), float3(-0.7, 0, 0) },
        { float3(0, -0.7, 0), float3(0, 0.7, 0) },
        { float3(0, 0, 0), float3(0, 0, 0) },
        { float3(0, 0, 0.3), float3(0.6, 0.2, 0) }
    };
    // Metaball field radii of max influence
    float radii[N_METABALLS] = { 0.35, 0.35, 0.25, 0.35, 0.4 };
//
//#if N_METABALLS == 6
//    float3 keyFrameCenters[N_METABALLS][2] =
//    {
//        { float3(-0.7, 0, 0), float3(0.7,0, 0) },
//        { float3(0.7 , 0, 0), float3(-0.7, 0, 0) },
//        { float3(0, -0.7, 0), float3(0, 0.7, 0) },
//        { float3(0, 0.7, 0), float3(0, -0.7, 0) },
//        { float3(0, 0, 0), float3(0, 0, 0) },
//        { float3(0, 0, 0.3), float3(0.6, 0.2, 0) }
//    };
//    // Metaball field radii of max influence
//    float radii[N_METABALLS] = { 0.35, 0.35, 0.35, 0.35, 0.25, 0.4 };

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
void FindIntersectingMetaballs(in Ray ray, out float tmin, out float tmax, inout Metaball blobs[N_METABALLS], \
    inout float centerDistance[N_METABALLS], out float blobsTmin[N_METABALLS], out float blobsTmax[N_METABALLS], \
    out UINT nActiveMetaballs)
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
            blobsTmin[nActiveMetaballs] = _thit;
            blobsTmax[nActiveMetaballs] = _tmax;
            nActiveMetaballs++;
#else
            nActiveMetaballs = N_METABALLS;
            blobsTmin[i] = _thit;
            blobsTmax[i] = _tmax;
#endif
        }
    }
    tmin = max(tmin, RayTMin());
    tmax = min(tmax, RayTCurrent());
}

// This function uses STATIC LOOPS, is called before FindIntersectingMetaballs().
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

void DivideGroups(in float centerDistance[N_METABALLS], in Metaball blobs[N_METABALLS], \
    in float blobsTmin[N_METABALLS], in float blobsTmax[N_METABALLS], \
    out Group groups[N_GROUPS], in UINT nActiveMetaballs, out UINT groupNum)
{
    // Number of non-empty groups
    groupNum = 0;

#if USE_DYNAMIC_LOOPS
    for (UINT i = 0; i < nActiveMetaballs; i++)
#else
    for (UINT i = 0; i < N_METABALLS; i++)
#endif
    {
        // Whether the new blob should be added to current group
        if ((groupNum > 0) && (((centerDistance[i] - blobs[i].radius) < groups[groupNum - 1].far) || (groupNum == N_GROUPS)))
        {
            // Append the new blob into current group
            groups[groupNum - 1].tail = i;
            groups[groupNum - 1].far = max(groups[groupNum - 1].far, (centerDistance[i] + blobs[i].radius));
            groups[groupNum - 1].tmin = min(groups[groupNum - 1].tmin, blobsTmin[i]);
            groups[groupNum - 1].tmax = max(groups[groupNum - 1].tmax, blobsTmax[i]);

            if ((centerDistance[i] - blobs[i].radius) < groups[groupNum - 1].near)
            { 
                // The new blob refreshes the group's `near`
                groups[groupNum - 1].near = centerDistance[i] - blobs[i].radius;

                while (groupNum > 1)
                { 
                    // Whether the new `near` leads to overlapping of former groups
                    if (groups[groupNum - 2].far > groups[groupNum - 1].near)
                    { 
                        // Decrement group counter
                        --groupNum;
                        // Merge two groups
                        groups[groupNum - 1].tail = groups[groupNum].tail;
                        groups[groupNum - 1].far = groups[groupNum].far;
                        groups[groupNum - 1].tmin = min(groups[groupNum - 1].tmin, groups[groupNum].tmin);
                        groups[groupNum - 1].tmax = max(groups[groupNum - 1].tmax, groups[groupNum].tmax);
                    }
                }
            }
        }

        // Create a new group
        else
        {
            groups[groupNum].head = i;
            groups[groupNum].tail = i;
            groups[groupNum].near = centerDistance[i] - blobs[i].radius;
            groups[groupNum].far = centerDistance[i] + blobs[i].radius;
            groups[groupNum].tmin = blobsTmin[i];
            groups[groupNum].tmax = blobsTmax[i];
            ++groupNum;
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
    
    float blobsTmin[N_METABALLS];
    float blobsTmax[N_METABALLS];
    FindIntersectingMetaballs(ray, tmin, tmax, blobs, centerDistance, blobsTmin, blobsTmax, nActiveMetaballs);


    Group groups[N_GROUPS];
    UINT groupNum = 0;
    DivideGroups(centerDistance, blobs, blobsTmin, blobsTmax, groups, nActiveMetaballs, groupNum);

    //**test**by zsym
    UINT MAX_STEPS = 64;

    for (UINT curGroupIdx = 0; curGroupIdx < groupNum; ++curGroupIdx) {
        float groupTmin = groups[curGroupIdx].tmin;
        float groupTmax = groups[curGroupIdx].tmax;
        float t = groupTmin;
        float minTStep = (groupTmax - groupTmin) / (MAX_STEPS / 1);
        UINT iStep = 0;

        while (iStep++ < MAX_STEPS)
        {
            float3 position = ray.origin + t * ray.direction;
            float fieldPotentials[N_METABALLS];    // Field potentials for each metaball.
            float sumFieldPotential = 0;           // Sum of all metaball field potentials.

            // Calculate field potentials from all metaballs in current group
            for (UINT i = groups[curGroupIdx].head; i <= groups[curGroupIdx].tail; ++i)
            {
                float distance;
                fieldPotentials[i] = CalculateMetaballPotential(position, blobs[i], distance);
                sumFieldPotential += fieldPotentials[i];
            }

            // Field potential threshold defining the isosurface.
            // Threshold - valid range is (0, 1), the larger the threshold the smaller the blob.
            const float Threshold = 0.25f;

            // Have we crossed the isosurface?
            if (sumFieldPotential >= Threshold)
            {
                //float3 normal = CalculateMetaballsNormal(position, blobs, nActiveMetaballs);
                float3 normal = CalculateGroupMetaballsNormal(position, blobs, groups[curGroupIdx].head, groups[curGroupIdx].tail);
                if (IsAValidHit(ray, t, normal))
                {
                    thit = t;
                    attr.normal = normal;
                    return true;
                }
            }

            t += minTStep;
        }
    }

    return false;
}

#endif // VOLUMETRICPRIMITIVESLIBRARY_H