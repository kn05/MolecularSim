# include <Siv3D.hpp> // OpenSiv3D v0.4.2
#include <Windows.h>

double dt=0.1;
double mass = 10;
double G = 1;

class vec2 {
private:
    double x, y;

public:
    vec2() :x(0.0), y(0.0) {}
    vec2(double x_, double y_) {
        x = x_;
        y = y_;
    }

    double getx() {
        return x;
    }
    double gety() {
        return y;
    }

    vec2 operator + (vec2 p) {
        double x_ = x + p.x;
        double y_ = y + p.y;
        return vec2(x_, y_);
    }
    vec2 operator - (vec2 p) {
        double x_ = x - p.x;
        double y_ = y - p.y;
        return vec2(x_, y_);
    }
    vec2 operator * (double p) {
        double x_ = x * p;
        double y_ = y * p;
        return vec2(x_, y_);
    }
    vec2 operator / (double p) {
        double x_ = x / p;
        double y_ = y / p;
        return vec2(x_, y_);
    }

    vec2 operator = (const vec2& p) {
        x = p.x;
        y = p.y;
        return vec2(x, y);
    }

    vec2 setlength(double p) {
        double d2 = pow(x, 2) + pow(y, 2);
        double d = pow(d2, 0.5);
        x = x * p / d;
        y = y * p / d;
        return vec2(x, y);
    }

    double length() {
        double d2 = pow(x, 2) + pow(y, 2);
        double d = pow(d2, 0.5);
        return d;
    }
};

class Molecular {
private:

public:
	Molecular() {

	};
	void draw() {

	};
	Vec2 pos;
	Vec2 vel;
};

using PVector = std::pair<std::vector<vec2>, std::vector<vec2>>;

std::pair<std::vector<vec2>, std::vector<vec2> > rk4(std::vector<vec2>& pos, std::vector<vec2>& vel);
std::vector<vec2> cal_vec(std::vector<vec2>& v0, std::vector<vec2>& k, int a);
std::vector<vec2> cal_vec(std::vector<vec2>& v0, std::vector<std::vector<vec2>>& k);
std::vector<vec2> cal_acc(std::vector<vec2>& pos);

void Main()
{
    std::vector<vec2> pos;
    std::vector<vec2> vel;
    pos.push_back(vec2(300, 400));
    pos.push_back(vec2(400, 400));
    vel.push_back(vec2(0, sqrt(2 * G * mass / 100)));
    vel.push_back(vec2(0, -sqrt(2 * G * mass / 100)));
	while (System::Update())
	{
        PVector t = std::make_pair(pos, vel);
		t = rk4(pos, vel);
        pos = t.first;
        vel = t.second;
        Print<< pos[1].getx() <<U" "<<pos[1].gety() ;
	}

}

PVector rk4(std::vector<vec2>& pos, std::vector<vec2>& vel) {
	std::vector<std::vector<vec2>> k1(4), k2(4);
    std::vector<std::vector<vec2>> x(4), v(4);
	x[0] = pos;
	v[0] = vel;
	k1[0] = v[0];
	k2[0] = cal_acc(x[0]);
    x[1] = cal_vec(x[0], k1[0], 0);
    v[1] = cal_vec(v[0], k2[0], 0);
    k1[1] = v[1];
    k2[1] = cal_acc(x[1]);
    x[2] = cal_vec(x[0], k1[1], 0);
    v[2] = cal_vec(v[0], k2[1], 0);
    k1[2] = v[2];
    k2[2] = cal_acc(x[2]);
    x[2] = cal_vec(x[0], k1[2], 1);
    v[2] = cal_vec(v[0], k2[2], 1);
    k1[3] = v[3];
    k2[3] = cal_acc(x[3]);
    x[3] = cal_vec(x[0], k1);
    v[3] = cal_vec(v[0], k2);
    return std::make_pair(x[3], v[3]);
}

std::vector<vec2> cal_vec(std::vector<vec2>& v0, std::vector<vec2>& k, int a) {
    int n = v0.size();
    std::vector<vec2> v(n);
    for (int i = 0; i < n; i++) {
        v[i] = v0[i] + k[i] * dt;
        if (!a) {
            v[i] = v[i] / 2;
        }
    }
    return v;
}
std::vector<vec2> cal_vec(std::vector<vec2>& v0, std::vector<std::vector<vec2>>& k) {
    int n = v0.size();
    std::vector<vec2> v(n);
    for (int i = 0; i < n; i++) {
        v[i] = v0[i] + (k[0][i] + k[1][i] * 2 + k[2][i] * 2 + k[0][i]) * dt / 6;
    }
    return v;
}

std::vector<vec2> cal_acc(std::vector<vec2>& pos) {
    int n = pos.size();
    std::vector<vec2> acc;
    for (auto &i:pos) {
        vec2 force = vec2(0, 0);
        for (auto &j:pos) {
            vec2 g = j - i;
            double d = g.length();
            if (d != 0) {
                g.setlength(G * mass * mass / (d * d));
                force = force + g;
            }
        }
        acc.push_back(force / mass);
    }
    return acc;
}
