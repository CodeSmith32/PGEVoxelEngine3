#ifndef OLC_PGEX_VOXELENGINE_H
#define OLC_PGEX_VOXELENGINE_H

#include "olcPixelGameEngine.h"

#include <unordered_map>

#define WORLD_ZMAX 256

namespace olc {

    namespace vox {

        struct Voxel {
            uint32_t colors[WORLD_ZMAX];
        };

        struct State {
            // world
            Voxel *world;
            int worldWidth;
            int worldHeight;
            bool worldOptimized;
            uint32_t background;

            // Screen
            int screenWidth;
            int screenHeight;
            int halfX;
            int halfY;
            float *depthBuffer;
            int16_t *jumpBuffer;
            float aspect;

            // Camera
            float fov;
            float fovXLength;
            float fovYLength;
            float renderDistance;
            float mipDistance;
            float cameraX;
            float cameraY;
            float cameraZ;
            float cameraDir;
            float cameraUp;
            float camCos;
            float camSin;
            float camTanUp;
        };

        struct ObjectView {
            float vx;
            float vy;
            float vz;
            float camDist;
            float camOff;
            float camAngle;
            float camUpAngle;
            int screenX;
            int screenY;
            int screenW;
            int screenH;
        };

        struct HitDetails {
            // if hit
            bool hit;
            // clamp percent of distance along vector
            float T;
            // new point
            float x;
            float y;
            float z;
            // voxel hit normal
            float nx;
            float ny;
            float nz;
            // coordinate of hit voxel
            int hitX;
            int hitY;
            int hitZ;
        };

        class Object {
        public:
            float x, y, z;
            float width, height;

            Object(): x(0.f), y(0.f), z(0.f), width(1.f), height(1.f) {}
            Object(float x, float y, float z, float width, float height):
                x(x), y(y), z(z), width(width), height(height) {}
            virtual ~Object() {}

            virtual uint32_t* GetPixelBuffer(int &width, int &height, const ObjectView* details) = 0;
            virtual uint32_t AdjustPixel(uint32_t pixel, int x, int y, const ObjectView* details);
        };

        class StaticObject: public Object {
        public:
            StaticObject(): Object() {}
            StaticObject(float x, float y, float z, float width, float height):
                Object(x, y, z, width, height) {}
            virtual ~StaticObject() {}
        
            virtual uint32_t* GetPixelBuffer(int &width, int &height, const ObjectView* details);

            // set array of images:
            // - width, height: size of each image
            // - images: num of images
            // - imageArray: pointer to array (of length 'images') of pointers to pixel data arrays, one for each image
            void SetImages(int width, int height, int images, uint32_t **imageArray);
            float objectDir;

        private:
            int w;
            int h;
            int images;
            uint32_t **imageArray;
        };

        class Engine : public olc::PGEX {
        public:
            Engine(const int screen_w, const int screen_h, const float fov, const float renderDistance, const float mipDistance = 0);
            virtual ~Engine();

            void SetCamera(float x, float y, float z, float dir, float up);  // 'up' causes distortion; use in moderation
            float GetFOV() { return state->fov; }
            void SetFOV(float fov);
            void RayCast(HitDetails &hit, float x, float y, float z, float vx, float vy, float vz, float xs=0.f, float ys=0.f, float zs=0.f);
            void Render();

            std::unordered_map<uint32_t, Object*> mapObjects;

        protected:
            Voxel* InitWorldVoxels(int width, int height);
            void OptimizeWorldVoxels();
            virtual uint32_t AdjustVoxel(uint32_t color, const Voxel &voxel, int tile_x, int tile_y, int tile_z, float distance);
            void SetBackground(uint32_t back);
            void CastRay(int x);
            void RenderObject(Object* obj);
            void RasterPixel(int x, int y, float z, const uint32_t &color);

        private:
            ObjectView* objView;
            State* state;
        };
    }
}

#define EPS 1e-5
#define PI 3.14159265359f
#define PI2 6.28318530718f
#define PI12 1.570796326795f
#define DTR 0.01745329251994f
#define RTD 57.29577951308f
#define VOX_ISDRAWABLE(V) (((V) & 0x80000000) != 0)
#define VOX_ISCOLLIDABLE(V) (((V) & 0x40000000) != 0)

// #ifdef OLC_PGEX_VOXELENGINE
#undef OLC_PGEX_VOXELENGINE

uint32_t olc::vox::Object::AdjustPixel(uint32_t pixel, int x, int y, const ObjectView* details) {
    return pixel;
}

uint32_t* olc::vox::StaticObject::GetPixelBuffer(int &width, int &height, const ObjectView* details) {
    width = w;
    height = h;
    float trueAngle = atan2f(-details->vy, -details->vx) - objectDir;

    int img = int((trueAngle + PI2) / PI2 * float(images) + 0.5f) % images;

    return imageArray[img];
}

void olc::vox::StaticObject::SetImages(int width, int height, int images, uint32_t **imageArray) {
    w = width;
    h = height;
    this->images = images;
    this->imageArray = imageArray;
}

olc::vox::Engine::Engine(const int screen_w, const int screen_h, const float fov, const float renderDistance, const float mipDistance) {
    objView = new olc::vox::ObjectView;
    state = new olc::vox::State;

    // init world state
    state->world = nullptr;
    state->worldOptimized = false;
    state->background = 0ul;

    // init screen state
    state->screenWidth = screen_w;
    state->screenHeight = screen_h;
    state->halfX = state->screenWidth / 2;
    state->halfY = state->screenHeight / 2;
    state->depthBuffer = new float[screen_w * screen_h];
    state->jumpBuffer = new int16_t[screen_w * screen_h];
    state->aspect = float(state->screenWidth) / float(state->screenHeight);

    // init camera state
    SetFOV(fov);
    state->renderDistance = renderDistance;
    state->mipDistance = mipDistance;
}

olc::vox::Engine::~Engine() {
    delete objView;

    if(state->world != nullptr)
        delete[] state->world;

    delete[] state->depthBuffer;

    delete state;
}

olc::vox::Voxel* olc::vox::Engine::InitWorldVoxels(int width, int height) {
    state->worldWidth = width;
    state->worldHeight = height;
    state->world = new olc::vox::Voxel[width * height];

    return state->world;
}

void olc::vox::Engine::OptimizeWorldVoxels() {
    for(int y=0; y<state->worldHeight; y++)
    for(int x=0; x<state->worldWidth; x++) {
        olc::vox::Voxel *vox = state->world + x + y * state->worldWidth;

        int solidZ = -1;
        for(int z=0; z<WORLD_ZMAX; z++) {
            if(VOX_ISDRAWABLE(vox->colors[z]))
                solidZ = z;
            else
                vox->colors[z] = (vox->colors[z] & 0xff000000) | (uint32_t(z - solidZ) & 0x00ffffff);
        }
    }
    state->worldOptimized = true;
}

uint32_t olc::vox::Engine::AdjustVoxel(uint32_t color, const olc::vox::Voxel& voxel, int tile_x, int tile_y, int tile_z, float distance) {
    return color;
}

void olc::vox::Engine::SetBackground(uint32_t back) {
    state->background = back;
}

void olc::vox::Engine::RayCast(olc::vox::HitDetails &hit, float x, float y, float z, float vx, float vy, float vz, float xs, float ys, float zs) {
    xs /= 2.f;
    ys /= 2.f;
    zs /= 2.f;

    // size cannot be absolute 0 since 1d rays can fit between voxels
    if(xs < EPS * 10) xs = EPS * 10;
    if(ys < EPS * 10) ys = EPS * 10;
    if(zs < EPS * 10) zs = EPS * 10;

    // avoid issues when vectors are too close but not equal to 0
    if(fabs(vx) < EPS) vx = 0;
    if(fabs(vy) < EPS) vy = 0;
    if(fabs(vz) < EPS) vz = 0;

    // vector heading step
    int stepX = vx > 0.f ? 1 : -1;
    int stepY = vy > 0.f ? 1 : -1;
    int stepZ = vz > 0.f ? 1 : -1;

    // offset to head-on corner
    float outX = vx >= 0.f ? xs : -xs;
    float outY = vy >= 0.f ? ys : -ys;
    float outZ = vz >= 0.f ? zs : -zs;

    // head-on corner position
    float px = x + outX;
    float py = y + outY;
    float pz = z + outZ;

    // cell position at head-on corner
    float nextX = floorf(px);
    float nextY = floorf(py);
    float nextZ = floorf(pz);

    // cell position
    int tileX = int(nextX);
    int tileY = int(nextY);
    int tileZ = int(nextZ);

    int xFromSide = vx == 0.f && nextX == px;
    int yFromSide = vy == 0.f && nextY == py;
    int zFromSide = vz == 0.f && nextZ == pz;

    // for edge case, set cell to current in order to test the immediate surface
    if(vx >= 0.f && nextX == px) { nextX--; tileX--; }
    if(vy >= 0.f && nextY == py) { nextY--; tileY--; }
    if(vz >= 0.f && nextZ == pz) { nextZ--; tileZ--; }

    // move nextXYZ to head-on corner of cell from corner of block
    if(vx > 0.f) nextX++;
    if(vy > 0.f) nextY++;
    if(vz > 0.f) nextZ++;

    float T = 0;
    Voxel* vox;

    // iterate voxel space
    while(T < 1) {
        // compute ray vector steps till next x / y / z plane
        float xf = fabs(vx) < EPS ? INFINITY : (nextX - px) / vx;
        float yf = fabs(vy) < EPS ? INFINITY : (nextY - py) / vy;
        float zf = fabs(vz) < EPS ? INFINITY : (nextZ - pz) / vz;

        if(xf < yf && xf < zf) { // x-hit is closest
            if(T + xf > 1) break;

            px = nextX;
            py += vy * xf;
            pz += vz * xf;

            nextX += stepX;
            tileX += stepX;
            T += xf;

            if(tileX < 0 || tileX >= state->worldWidth) continue;

            // test y-z surface for non-empty voxels
            int y1 = int(floorf(py - outY - ys));
            int y2 = int(floorf(py - outY + ys)) + 1 - yFromSide;
            int z1 = int(floorf(pz - outZ - zs));
            int z2 = int(floorf(pz - outZ + zs)) + 1 - zFromSide;

            if(y1 < 0) y1 = 0; else if(y1 > state->worldHeight) y1 = state->worldHeight;
            if(y2 < 0) y2 = 0; else if(y2 > state->worldHeight) y2 = state->worldHeight;
            if(z1 < 0) z1 = 0; else if(z1 > WORLD_ZMAX) z1 = WORLD_ZMAX;
            if(z2 < 0) z2 = 0; else if(z2 > WORLD_ZMAX) z2 = WORLD_ZMAX;

            for(int y=y1; y<y2; y++) {
                vox = state->world + tileX + y*state->worldWidth;
                
                for(int z=z1; z<z2; z++) {
                    if(VOX_ISCOLLIDABLE(vox->colors[z])) {
                        hit.hit = true;
                        hit.T = T;
                        hit.x = px - outX;
                        hit.y = py - outY;
                        hit.z = pz - outZ;
                        hit.nx = -stepX;
                        hit.ny = 0;
                        hit.nz = 0;
                        hit.hitX = tileX;
                        hit.hitY = y;
                        hit.hitZ = z;
                        return;
                    }
                }
            }
        } else if(yf < zf) { // y-hit is closest
            if(T + yf > 1) break;

            px += vx * yf;
            py = nextY;
            pz += vz * yf;

            nextY += stepY;
            tileY += stepY;
            T += yf;

            if(tileY < 0 || tileY >= state->worldWidth) continue;

            // test x-z surface for non-empty voxels
            int x1 = int(floorf(px - outX - xs));
            int x2 = int(floorf(px - outX + xs)) + 1 - xFromSide;
            int z1 = int(floorf(pz - outZ - zs));
            int z2 = int(floorf(pz - outZ + zs)) + 1 - zFromSide;

            if(x1 < 0) x1 = 0; else if(x1 > state->worldWidth) x1 = state->worldWidth;
            if(x2 < 0) x2 = 0; else if(x2 > state->worldWidth) x2 = state->worldWidth;
            if(z1 < 0) z1 = 0; else if(z1 > WORLD_ZMAX) z1 = WORLD_ZMAX;
            if(z2 < 0) z2 = 0; else if(z2 > WORLD_ZMAX) z2 = WORLD_ZMAX;

            for(int x=x1; x<x2; x++) {
                vox = state->world + x + tileY*state->worldWidth;

                for(int z=z1; z<z2; z++) {
                    if(VOX_ISCOLLIDABLE(vox->colors[z])) {
                        hit.hit = true;
                        hit.T = T;
                        hit.x = px - outX;
                        hit.y = py - outY;
                        hit.z = pz - outZ;
                        hit.nx = 0;
                        hit.ny = -stepY;
                        hit.nz = 0;
                        hit.hitX = tileX;
                        hit.hitY = y;
                        hit.hitZ = z;
                        return;
                    }
                }
            }
        } else {
            if(T + zf > 1) break;

            px += vx * zf;
            py += vy * zf;
            pz = nextZ;

            nextZ += stepZ;
            tileZ += stepZ;
            T += zf;

            if(tileZ < 0 || tileZ >= WORLD_ZMAX) continue;

            // test x-y surface for non-empty voxels
            int x1 = int(floorf(px - outX - xs));
            int x2 = int(floorf(px - outX + xs)) + 1 - xFromSide;
            int y1 = int(floorf(py - outY - ys));
            int y2 = int(floorf(py - outY + ys)) + 1 - yFromSide;
            
            if(x1 < 0) x1 = 0; else if(x1 > state->worldWidth) x1 = state->worldWidth;
            if(x2 < 0) x2 = 0; else if(x2 > state->worldWidth) x2 = state->worldWidth;
            if(y1 < 0) y1 = 0; else if(y1 > state->worldHeight) y1 = state->worldHeight;
            if(y2 < 0) y2 = 0; else if(y2 > state->worldHeight) y2 = state->worldHeight;

            for(int y=y1; y<y2; y++)
            for(int x=x1; x<x2; x++) {
                vox = state->world + x + y*state->worldWidth;

                if(VOX_ISCOLLIDABLE(vox->colors[tileZ])) {
                    hit.hit = true;
                    hit.T = T;
                    hit.x = px - outX;
                    hit.y = py - outY;
                    hit.z = pz - outZ;
                    hit.nx = 0;
                    hit.ny = 0;
                    hit.nz = -stepZ;
                    hit.hitX = x;
                    hit.hitY = y;
                    hit.hitZ = tileZ;
                    return;
                }
            }
        }
    }

    // no hit
    hit.hit = false;
    hit.T = 1;
    hit.x = x + vx;
    hit.y = y + vy;
    hit.z = z + vz;
    hit.nx = 0;
    hit.ny = 0;
    hit.nz = 0;
    hit.hitX = -1;
    hit.hitY = -1;
    hit.hitZ = -1;
}

void olc::vox::Engine::SetCamera(float x, float y, float z, float dir, float up) {
    state->cameraX = x;
    state->cameraY = y;
    state->cameraZ = z;
    state->cameraDir = dir;
    state->cameraUp = up; // causes distortion; use in moderation
}

void olc::vox::Engine::SetFOV(float fov) {
    state->fov = fov;
    state->fovYLength = tanf(fov * 0.5f);
    state->fovXLength = state->fovYLength * state->aspect;
}

void olc::vox::Engine::Render() {
    state->camCos = cos(state->cameraDir);
    state->camSin = sin(state->cameraDir);
    state->camTanUp = tan(state->cameraUp);

    if(!state->worldOptimized)
        OptimizeWorldVoxels();

    for(int i=state->screenWidth * state->screenHeight - 1; i>=0; i--) {
        state->depthBuffer[i] = INFINITY;
        state->jumpBuffer[i] = 0;
    }

    for(int x=0; x<state->screenWidth; x++) {
        CastRay(x);
    }

    objView->camAngle = state->cameraDir;
    objView->camUpAngle = state->cameraUp;

    for(const auto& objPair : mapObjects) {
        RenderObject(objPair.second);
    }
}

void olc::vox::Engine::CastRay(int x) {
    float fovOff = (1.f - float(x) / float(state->screenWidth) * 2.f) * state->fovXLength;
    int bufI = x * state->screenHeight; // column index

    // column cullers
    int p1 = state->screenHeight; // cull below p1
    int p3 = 0; // cull above p3
    int pa = 0, pb = 0; // cull between pa and pb
    bool afterAB = false; // if pa -> pb culling is past

    // ray tracing
    float d = 0; // distance traveled
    float u = sqrtf(fovOff * fovOff + 1.f); // length of ray vector 'v'
    float px = state->cameraX;
    float py = state->cameraY; // current ray test point
    float vx = state->camCos - fovOff*state->camSin;
    float vy = state->camSin + fovOff*state->camCos; // ray vector (projected to rotated surface)
    int stepX = vx > 0 ? 1 : -1;
    int stepY = vy > 0 ? 1 : -1; // step direction of ray vector
    float nextX = float(floor(px));
    float nextY = float(floor(py)); // point on 'tile' that ray vector is facing
    int tileX = int(nextX);
    int tileY = int(nextY); // the voxel / tile that the ray is currently in

    // jump point to the corner of 'tile' that the ray is facing
    if(vx > 0.f) nextX++;
    if(vy > 0.f) nextY++;

    // edge case: 
    if(vx < 0.f && nextX == px) { nextX--; tileX--; }
    if(vy < 0.f && nextY == py) { nextY--; tileY--; }

    float nextDistance = state->mipDistance;
    float halfCell = 0.5f;
    int cellSize = 1;

    while(d < state->renderDistance) {
        if(d > nextDistance && cellSize < WORLD_ZMAX) {
            nextDistance *= 2;

            stepX *= 2; // make step 2x2
            stepY *= 2; // make step 2x2
            halfCell *= 2; // make half-cell for center of 2x2
            cellSize *= 2; // double cellsize
            tileX = tileX / cellSize * cellSize; // truncate tile to 2x2
            tileY = tileY / cellSize * cellSize; // truncate tile to 2x2

            // update next
            nextX = tileX;
            nextY = tileY;
            if(vx > 0.f) nextX += cellSize;
            if(vy > 0.f) nextY += cellSize;

            // edge case:
            if(vx < 0.f && nextX == px) { nextX -= cellSize; tileX -= cellSize; }
            if(vy < 0.f && nextY == py) { nextY -= cellSize; tileY -= cellSize; }
        }

        // compute ray vector steps till next x / y plane
        float xf = fabs(vx) < EPS ? INFINITY : (nextX - px) / vx;
        float yf = fabs(vy) < EPS ? INFINITY : (nextY - py) / vy;

        if(xf < 0.f || yf < 0.f) break; // sanity check; shouldn't ever happen

        if(xf < yf) { // x-intersection is closer
            // step cell
            nextX += stepX;
            tileX += stepX;

            d += u * xf; // increase distance
            px += vx * xf; // progress ray
            py += vy * xf; // progress ray

            if(stepX < 0 ? tileX < 0 : tileX >= state->worldWidth) break; // ray is headed outside room
        } else { // y-intersection is closer
            // step cell
            nextY += stepY;
            tileY += stepY;

            d += u * yf; // increase distance
            px += vx * yf; // progress ray
            py += vy * yf; // progress ray

            if(stepY < 0 ? tileY < 0 : tileY >= state->worldHeight) break; // ray is headed outside room
        }

        if(tileX < 0 || tileX >= state->worldWidth || tileY < 0 || tileY >= state->worldHeight) continue;

        olc::vox::Voxel* vox = state->world + tileX + tileY * state->worldWidth;

        float dd = ( // dot-product projection distance
            (float(tileX) + halfCell - state->cameraX) * state->camCos
            + (float(tileY) + halfCell - state->cameraY) * state->camSin
        );
        if(dd <= EPS) continue;

        // column render process

        float projZ = state->cameraZ + dd * state->camTanUp;
        float halfZ = dd * state->fovYLength;

        // cast from screen to world
        #define SCREEN_TO_WORLD(Y) \
            int(projZ + (1.0f - float(Y) / float(state->halfY)) * halfZ)

        // cast from world to screen
        #define WORLD_TO_SCREEN(Z) \
            int((1.0f - (float(Z) - projZ) / halfZ) * state->halfY)

        int y = 0, jumpFromY, nextY, prevY, minY, jumpY;
        float voxScale;
        bool drawn = false, setjump = false;
        uint32_t color;

        int tileZ = SCREEN_TO_WORLD(0.0f),
            endZ = SCREEN_TO_WORLD(state->screenHeight-1) / cellSize * cellSize;
        if(tileZ < 0 || endZ >= WORLD_ZMAX) goto BAIL; // bail: out of view cuz we're too high up or too low underground
        tileZ = tileZ / cellSize * cellSize; // truncate tileZ after < 0 check since truncation might truncate it 'up' to 0

        if(tileZ >= WORLD_ZMAX) {
            y = WORLD_TO_SCREEN(WORLD_ZMAX);
            if(y < 0) y = 0;
            if(y >= state->screenHeight) goto BAIL;
            tileZ = WORLD_ZMAX - cellSize;
        }
        if(endZ < 0) endZ = 0;

        // steep angle down voxel-mortar: if offscreen voxel at top scales into view, then start with it
        if(tileZ + cellSize < WORLD_ZMAX && VOX_ISDRAWABLE(vox->colors[tileZ + cellSize])) {
            voxScale = ((state->cameraZ - float(tileZ + cellSize)) / dd - 1.0f) * float(cellSize) * 2.0f;
            if(voxScale > 0.0f && WORLD_TO_SCREEN(float(tileZ + cellSize) - voxScale) + 1 >= 0)
                tileZ += cellSize;
        }
        // steep angle up voxel-mortar: if offscreen voxel at bottom scales into view, then include it at end
        if(endZ - cellSize >= 0 && VOX_ISDRAWABLE(vox->colors[endZ - cellSize])) {
            voxScale = ((float(endZ - cellSize) - state->cameraZ) / dd - 1.0f) * float(cellSize) * 2.0f;
            if(voxScale > 0.0f && WORLD_TO_SCREEN(float(endZ - cellSize) + voxScale) < state->screenHeight)
                endZ -= cellSize;
        }

        // render column of screen pixels
        while(1) {
            // handle empty voxel jumps
            if(!VOX_ISDRAWABLE(vox->colors[tileZ])) {
                do {
                    tileZ -= int(vox->colors[tileZ]);
                    if(tileZ < endZ) goto BAIL; // we've passed the end voxel; nothing else opaque to render
                    tileZ = tileZ / cellSize * cellSize;
                } while(!VOX_ISDRAWABLE(vox->colors[tileZ]));

                prevY = y + 1;

                // steep angle up voxel-mortar
                voxScale = fmax((float(tileZ + cellSize) - state->cameraZ) / dd - 1.0f, 0.0f) * float(cellSize) * 2.0f;
                y = WORLD_TO_SCREEN((tileZ + cellSize) + voxScale); // '+ voxScale' fills in steep projection angle gaps
                if(y < prevY) y = prevY;
            }

            // handle buffer jumps
            if(state->jumpBuffer[bufI+y] != 0) {
                minY = y;
                while((jumpY = state->jumpBuffer[bufI+y]) != 0) {
                    y += jumpY; // keep jumping
                    if(y < minY) minY = y;

                    if(y >= state->screenHeight) {
                        state->jumpBuffer[bufI+minY] = state->screenHeight - minY;
                        goto BAIL;
                    }
                }
                state->jumpBuffer[bufI+minY] = y - minY;

                int prevZ = tileZ;
                tileZ = SCREEN_TO_WORLD(y);
                if(tileZ < endZ) goto BAIL;
                tileZ = tileZ / cellSize * cellSize;
                if(tileZ > prevZ) tileZ = prevZ;

                continue;
            }

            // steep angle down voxel-mortar
            voxScale = fmax((state->cameraZ - float(tileZ)) / dd - 1.0f, 0.0f) * float(cellSize) * 2.0f;
            nextY = WORLD_TO_SCREEN(float(tileZ) - voxScale); // '- voxScale' fills in steep projection angle gaps
            if(nextY < y) { // if voxel is too small and doesn't change y, keep going
                tileZ--;
                if(tileZ < endZ) goto BAIL;
                continue;
            }

            // render single voxel column
            color = AdjustVoxel(vox->colors[tileZ], *vox, tileX, tileY, tileZ, d);
            jumpFromY = y;
            while(1) {
                RasterPixel(x, y, dd, color);
                state->jumpBuffer[bufI+y] = jumpFromY - y; // set pointbacks

                if(++y >= state->screenHeight) {
                    state->jumpBuffer[bufI+jumpFromY] = state->screenHeight - jumpFromY;
                    goto BAIL;
                }
                if(state->jumpBuffer[bufI+y] != 0) break; // hit a jump buffer
                if(y > nextY) { // finished voxel
                    tileZ--;
                    break;
                }
            }
            state->jumpBuffer[bufI+jumpFromY] = y - jumpFromY;
            if(tileZ < endZ) break;
        }

        BAIL:
        if(state->jumpBuffer[bufI] >= state->screenHeight) break;
    }

    // fill in sky
    int jumpY;
    for(int y=0; y<state->screenHeight; y++) {
        while((jumpY = state->jumpBuffer[bufI+y]) != 0) {
            y += jumpY;
            if(y >= state->screenHeight) goto END;
        }

        RasterPixel(x, y, INFINITY, state->background);
    }

    END:
    return;
}

void olc::vox::Engine::RenderObject(olc::vox::Object* obj) {
    // get cam to object vector
    objView->vx = obj->x - state->cameraX;
    objView->vy = obj->y - state->cameraY;
    objView->vz = obj->z - state->cameraZ;

    objView->camDist = objView->vx * state->camCos + objView->vy * state->camSin; // dot product
    if(objView->camDist < 0.5f) return; // behind or too close

    objView->camOff = objView->vx * state->camSin - objView->vy * state->camCos; // perp dot product (90deg CW from camera)
    float layerXSize = state->fovXLength * objView->camDist;
    float layerYSize = state->fovYLength * objView->camDist;
    float zOff = objView->vz - objView->camDist * state->camTanUp;
    if(fabs(objView->camOff) - layerXSize > obj->width * 0.5f || fabs(zOff) - layerYSize > obj->height * 0.5f) return; // outside camera fov

    float xFactor = float(state->screenWidth) * 0.5f / layerXSize;
    float zFactor = float(state->screenHeight) * 0.5f / layerYSize;

    objView->screenX = state->halfX + int(objView->camOff * xFactor);
    objView->screenY = state->halfY - int(zOff * zFactor);
    objView->screenW = int(obj->width * xFactor);
    objView->screenH = int(obj->height * zFactor);

    objView->screenX -= objView->screenW / 2;
    objView->screenY -= objView->screenH / 2;

    int width = 0, height = 0;
    uint32_t* pixels = obj->GetPixelBuffer(width, height, objView);
    if(pixels == nullptr || !width || !height) return;

    int x1 = objView->screenX < 0 ? 0 : objView->screenX;
    int y1 = objView->screenY < 0 ? 0 : objView->screenY;
    int x2 = objView->screenX + objView->screenW > state->screenWidth ? state->screenWidth : objView->screenX + objView->screenW;
    int y2 = objView->screenY + objView->screenH > state->screenHeight ? state->screenHeight : objView->screenY + objView->screenH;

    for(int y=y1; y<y2; y++)
    for(int x=x1; x<x2; x++) {
        int tx = int((float(x) - float(objView->screenX)) / float(objView->screenW) * float(width));
        int ty = int((float(y) - float(objView->screenY)) / float(objView->screenH) * float(height));

        uint32_t color = obj->AdjustPixel(pixels[ty*width + tx], tx, ty, objView);

        RasterPixel(x, y, objView->camDist, color);
    }
}

void olc::vox::Engine::RasterPixel(int x, int y, float z, const uint32_t &color) {
    int i = x*state->screenHeight + y;

    if(z <= state->depthBuffer[i]) {
        pge->Draw(x, y, color | 0xff000000);
        state->depthBuffer[i] = z;
    }
}

// #endif // OLC_PGEX_VOXELENGINE

#endif // OLC_PGEX_VOXELENGINE_H