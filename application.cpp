#include "olcPixelGameEngine.h"
#define OLC_PGEX_VOXELENGINE
#include "olcPGEX_VoxelEngine.h"

#include "clock.h"

#include <unordered_map>
#include <string>

class VoxelTest : public olc::vox::Engine {
public:
	VoxelTest(const int screen_w, const int screen_h, const float fov);
	virtual ~VoxelTest();

	int width() { return worldWidth; }
	int height() { return worldHeight; }

protected:
	void GenerateWorldVoxels();
	uint32_t AdjustVoxel(uint32_t color, const olc::vox::Voxel& voxel, int tile_x, int tile_y, int tile_z, float distance) override;

	int worldWidth;
	int worldHeight;
	float fogNear;
	float fogFar;
	uint32_t fogColor;

	olc::Sprite* mapColor;
	olc::Sprite* mapHeight;
};

VoxelTest::VoxelTest(const int screen_w, const int screen_h, const float fov):
		olc::vox::Engine(screen_w, screen_h, fov, 500.f, 200.f /*1000.f, 150.f*/) {
	fogNear = 200.f;
	fogFar = 1500.f;
	fogColor = 0xFFFFA030;

	mapColor = new olc::Sprite("./map.color.png");
	mapHeight = new olc::Sprite("./map.height.png");

	worldWidth = mapColor->width;
	worldHeight = mapColor->height;

	SetBackground(fogColor);
	GenerateWorldVoxels();
}
VoxelTest::~VoxelTest() {}

void VoxelTest::GenerateWorldVoxels() {
	olc::vox::Voxel* grid = InitWorldVoxels(worldWidth, worldHeight);
	olc::Pixel color, height;

	for(int y=0; y<worldHeight; y++)
	for(int x=0; x<worldWidth; x++) {
		int i = y*worldWidth + x;
		color = mapColor->GetPixel(x, worldHeight - y - 1);

		int off = rand()%8 - 4,
			r = int(color.r) + off,
			g = int(color.g) + off,
			b = int(color.b) + off;

		color.r = uint8_t(r < 0 ? 0 : r > 255 ? 255 : r);
		color.g = uint8_t(g < 0 ? 0 : g > 255 ? 255 : g);
		color.b = uint8_t(b < 0 ? 0 : b > 255 ? 255 : b);

		height = mapHeight->GetPixel(x, worldHeight - y - 1);

		for(int j=0; j<256; j++)
			grid[i].colors[j] = /*j % 2 == 0 &&*/ (j <= height.r || (j > height.g && j <= height.b)) ? color.n : 0x00000000;
	}
	OptimizeWorldVoxels();
}

uint32_t VoxelTest::AdjustVoxel(uint32_t color, const olc::vox::Voxel& voxel, int tile_x, int tile_y, int tile_z, float distance) {
	// return voxel.color;

	// float T = (distance - fogNear) / (fogFar - fogNear);
	// T = T < 0.f ? 0.f : T > 1.f ? 1.f : T;
	// olc::Pixel *voxColor = (olc::Pixel*)&voxel.color;
	// return (*voxColor * (1 - T) + fogColor * T).n;

	float T = (distance - fogNear) / (fogFar - fogNear);
	T = T < 0.f ? 0.f : T > 1.f ? 1.f : T;
	float _T = 1 - T;

	float r1 = float(color & 255),
		  g1 = float((color >> 8) & 255),
		  b1 = float((color >> 16) & 255),
		  r2 = float(fogColor & 255),
		  g2 = float((fogColor >> 8) & 255),
		  b2 = float((fogColor >> 16) & 255);

	return 0xFF000000
		+ ((int(b1*_T + b2*T) & 255) << 16)
		+ ((int(g1*_T + g2*T) & 255) << 8)
		+ (int(r1*_T + r2*T) & 255);
}

class SpriteManager {
public:
	virtual ~SpriteManager();

	void AddSprite(const std::string name, olc::Sprite* spr);
	void DropSprite(const std::string name);
	olc::Sprite* GetSprite(const std::string name);

private:
	std::unordered_map<std::string, olc::Sprite*> spriteMap;
};

SpriteManager::~SpriteManager() {
	for(const auto& sprPair : spriteMap) {
		delete sprPair.second;
	}
}

void SpriteManager::AddSprite(const std::string name, olc::Sprite* spr) {
	spriteMap.insert_or_assign(name, spr);
}
void SpriteManager::DropSprite(const std::string name) {
	spriteMap.erase(name);
}
olc::Sprite* SpriteManager::GetSprite(const std::string name) {
	return spriteMap.at(name);
}

SpriteManager spriteManager;

class FPSLimiter {
	float fps;
	Clock timer;
	double lastTime, timeDiff;

public:
	FPSLimiter(): fps(INFINITY), lastTime(-1e100), timeDiff(0) {}
	FPSLimiter(float fps): fps(fps), lastTime(-1e100), timeDiff(1000.f / fps) {}
	virtual ~FPSLimiter() {}
	
	void SetFPS(float fps) { this->fps = fps; timeDiff = 1000.f / fps; }
	float GetFPS() { return fps; }
	void Wait() {
		while(timer.getMilliseconds() - lastTime < timeDiff);
		lastTime = timer.getMilliseconds();
	}
};

class TestObject: public olc::vox::StaticObject {
	static uint32_t **spriteArray;
	static int sprWidth;
	static int sprHeight;
public:
	TestObject(float x, float y, float z);
};

uint32_t **TestObject::spriteArray = nullptr;
int TestObject::sprWidth = 0;
int TestObject::sprHeight = 0;

TestObject::TestObject(float x, float y, float z):
		StaticObject(x, y, z, 8.f, 8.f) {
	if(spriteArray == nullptr) {
		olc::Sprite* spr = spriteManager.GetSprite("default");

		sprWidth = spr->width;
		sprHeight = spr->height;

		spriteArray = new uint32_t*;
		*spriteArray = (uint32_t*)spr->GetData();
	}
	SetImages(sprWidth, sprHeight, 1, spriteArray);
}

class Game : public olc::PixelGameEngine {
public:
	Game(uint32_t w, uint32_t h, uint32_t px = 1, uint32_t py = 1) {
		sAppName = "Voxel Engine Test";

		fpsLimiter.SetFPS(50.f);

		if(Construct(w, h, px, py)) {
			Start();
		}
	}

public:
	bool OnUserCreate() {
		fov = 90;

		world = new VoxelTest(ScreenWidth(), ScreenHeight(), float(fov) * DTR);

		olc::Sprite* defaultSpr = new olc::Sprite(16, 16);
		olc::Pixel* px = defaultSpr->GetData();
		for(int y=0; y<defaultSpr->height; y++)
		for(int x=0; x<defaultSpr->width; x++)
			px[x + y*defaultSpr->width] = 0xFF800080;

		spriteManager.AddSprite("default", defaultSpr);

		TestObject* obj1 = new TestObject(0.5f, 0.5f, 6.f);
		TestObject* obj2 = new TestObject(10.5f, 0.5f, 6.f);
		TestObject* obj3 = new TestObject(0.5f, 5.5f, 6.f);
		world->mapObjects.insert_or_assign(1, obj1);
		world->mapObjects.insert_or_assign(2, obj2);
		world->mapObjects.insert_or_assign(3, obj3);

		xcam = 0.0f; // float(world->width()) / 2.f;
		ycam = 0.0f; // float(world->height()) / 2.f;
		zcam = 20.f;
		dcam = 0.f;
		ucam = 0.f;
		speed = 1.4f;
		jumpSpeed = 3.f;
		gravity = 0.4f;
		sizeX = 1.8f;
		sizeY = 1.8f;
		sizeZ = 3.8f;
		camHeight = 3.2f;
		vx = 0.f;
		vy = 0.f;
		vz = 0.f;
		onground = false;
		flying = true;

		return true;
	}
	
	bool OnUserUpdate(float elapsedTime) {
		float inputFactor;

		inputFactor = GetKey(olc::LEFT).bHeld - GetKey(olc::RIGHT).bHeld;
		dcam += 0.05 * inputFactor * (GetKey(olc::SHIFT).bHeld ? 0.25f : 1.f);

		inputFactor = GetKey(olc::HOME).bHeld - GetKey(olc::END).bHeld;
		ucam += 0.05f * inputFactor * (GetKey(olc::SHIFT).bHeld ? 0.25f : 1.f);
		ucam = ucam < -PI * 0.4f ? -PI * 0.4f : ucam > PI * 0.4f ? PI * 0.4f : ucam;

		if(GetKey(olc::ENTER).bPressed) {
			if(GetKey(olc::CTRL).bHeld) {
				dcam = roundf(dcam / PI12 * 2.f) * PI12 * 0.5f;
				ucam = 0.f;
			} else
				flying = !flying;
		}


		if(dcam < 0) dcam += PI2;
		if(dcam > PI2) dcam -= PI2;

		float xface = cos(dcam) * speed, yface = sin(dcam) * speed;

		//inputFactor = GetKey(olc::PGUP).bHeld - GetKey(olc::PGDN).bHeld;
		//speed += inputFactor * 0.05;
		//if(inputFactor != 0) speed = speed < 0.1 ? 0.1 : speed > 50 ? 50 : speed;

		if(!flying) vz -= gravity; // gravity

		vx *= 0.95f;
		vy *= 0.95f;
		vz *= 0.95f;

		onground = false;

		olc::vox::HitDetails hit;
		float mx = vx, my = vy, mz = vz;
		for(int _=0; _<3; _++) {
			world->RayCast(hit, xcam, ycam, zcam, mx, my, mz, sizeX, sizeY, sizeZ);

			// remove portion of motion vector that progressed
			mx -= hit.x - xcam;
			my -= hit.y - ycam;
			mz -= hit.z - zcam;

			// progress camera to new position
			xcam = hit.x;
			ycam = hit.y;
			zcam = hit.z;

			if(hit.hit) {
				if(hit.nx) mx = vx = 0; // if x hit, stop moving vx
				if(hit.ny) my = vy = 0; // if y hit, stop moving vy
				if(hit.nz) mz = vz = 0; // if z hit, stop moving vz
				if(hit.nz > 0) onground = true; // if hit floor, we're onground

				if(!mx && !my && !mz) break; // if hit stopped our motion, stop testing
			} else break; // no hit, we moved full distance
		}

		// move forward / backward
		inputFactor = GetKey(olc::UP).bHeld - GetKey(olc::DOWN).bHeld;
		vx += xface * inputFactor;
		vy += yface * inputFactor;

		// strafe left / right
		inputFactor = GetKey(olc::Z).bHeld - GetKey(olc::X).bHeld;
		vx += -yface * inputFactor;
		vy += xface * inputFactor;

		vx *= 0.5f;
		vy *= 0.5f;

		if(flying) {
			vz += speed * (GetKey(olc::SPACE).bHeld - GetKey(olc::SHIFT).bHeld);
			vz *= 0.5f;
		} else if(onground) {
			// jump
			if(GetKey(olc::SPACE).bHeld) vz = jumpSpeed;
		}

		inputFactor = GetKey(olc::Q).bHeld - GetKey(olc::A).bHeld;
		fov += inputFactor;
		fov = fov < 1 ? 1 : fov > 160 ? 160 : fov;

		world->SetFOV(float(fov) * DTR);
		world->SetCamera(xcam, ycam, zcam - sizeZ*0.5f + camHeight, dcam, ucam);
		world->Render();

		DrawString({4,4}, "FPS "+std::to_string(int(1.0f / elapsedTime)));
		DrawString({ScreenWidth()/3, 4}, "FOV "+std::to_string(fov));
		DrawString({ScreenWidth()*2/3, 4}, "DIR "+std::to_string(int(dcam * RTD)));

		if(GetKey(olc::ESCAPE).bHeld) return false;

		fpsLimiter.Wait();

		return true;
	}
	
private:
	VoxelTest* world;

	FPSLimiter fpsLimiter;
	int fov; // integer in degrees
	float xcam, ycam, zcam, dcam, ucam, speed, jumpSpeed;
	float vx, vy, vz;
	float sizeX, sizeY, sizeZ, camHeight;
	float gravity;
	bool onground, flying;
};

int main(){
	Game app(320, 240, 3, 3);
	return 0;
}