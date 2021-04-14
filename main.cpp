#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <cassert>
using namespace std;

struct PRNG{
    std::mt19937 engine;
};

void initGenerator(PRNG& generator) {
    // Создаём псевдо-устройство для получения случайного зерна.
    std::random_device device;
    // Получаем случайное зерно последовательности
    generator.engine.seed(device());
}

// Генерирует целое число в диапазоне [minValue, maxValue)
unsigned random(PRNG& generator, unsigned minValue, unsigned maxValue) {
    // Проверяем корректность аргументов
    assert(minValue < maxValue);

    // Создаём распределение
    std::uniform_int_distribution<unsigned> distribution(minValue, maxValue);

    // Вычисляем псевдослучайное число: вызовем распределение как функцию,
    //  передав генератор произвольных целых чисел как аргумент.
    return distribution(generator.engine);
}

// Генерирует число с плавающей точкой в диапазоне [minValue, maxValue)
float getRandomFloat(PRNG& generator, float minValue, float maxValue) {
    // Проверяем корректность аргументов
    assert(minValue < maxValue);

    // Создаём распределение
    std::uniform_real_distribution<float> distribution(minValue, maxValue);

    // Вычисляем псевдослучайное число: вызовем распределение как функцию,
    //  передав генератор произвольных целых чисел как аргумент.
    return distribution(generator.engine);
}

class Particle {
public:
    double mass;
    double K;
    double U;

    double x;
    double y;
    double z;
    double x_0;
    double y_0;
    double z_0;
    double x_prev;
    double y_prev;
    double z_prev;
    double v_x;
    double v_y;
    double v_z;

    double F_x;
    double F_y;
    double F_z;
};

void Update_coords(Particle *p, int number_of_particles, int cage, double dt) {
    double temp;

    for (int i = 0; i < number_of_particles; i++) {
        temp = p[i].x;
        p[i].x = 2 * p[i].x - p[i].x_prev + p[i].F_x / p[i].mass * dt * dt;
        p[i].x_0 -= p[i].x;
        temp = p[i].x - temp;
        if (p[i].x >= cage) {
            p[i].x = fmod(p[i].x, cage);
        } else if (p[i].x < 0) {
            p[i].x = cage + fmod(p[i].x, cage);
        }
        p[i].x_prev = p[i].x - temp;
        p[i].x_0 += p[i].x;

        temp = p[i].y;
        p[i].y = 2 * p[i].y - p[i].y_prev + p[i].F_y / p[i].mass * dt * dt;
        p[i].y_0 -= p[i].y;
        temp = p[i].y - temp;
        if (p[i].y >= cage) {
            p[i].y = fmod(p[i].y, cage);
        } else if (p[i].y < 0) {
            p[i].y = cage + fmod(p[i].y, cage);
        }
        p[i].y_prev = p[i].y - temp;
        p[i].y_0 += p[i].y;

        temp = p[i].z;
        p[i].z = 2 * p[i].z - p[i].z_prev + p[i].F_z / p[i].mass * dt * dt;
        p[i].z_0 -= p[i].z;
        temp = p[i].z - temp;
        if (p[i].z >= cage) {
            p[i].z = fmod(p[i].z, cage);
        } else if (p[i].z < 0) {
            p[i].z = cage + fmod(p[i].z, cage);
        }
        p[i].z_prev = p[i].z - temp;
        p[i].z_0 += p[i].z;

        // Intermediate velocity
        p[i].v_x += 0.5 * p[i].F_x / p[i].mass * dt;
        p[i].v_y += 0.5 * p[i].F_y / p[i].mass * dt;
        p[i].v_z += 0.5 * p[i].F_z / p[i].mass * dt;
    }
}

void Update_velocity(Particle *p, int number_of_particles, double dt) {
    // New velocity
    for (int i = 0; i < number_of_particles; i++) {
        p[i].v_x += 0.5 * p[i].F_x / p[i].mass * dt;
        p[i].v_y += 0.5 * p[i].F_y / p[i].mass * dt;
        p[i].v_z += 0.5 * p[i].F_z / p[i].mass * dt;

        p[i].K = (p[i].v_x*p[i].v_x + p[i].v_y*p[i].v_y + p[i].v_z*p[i].v_z) / 2;
    }
}

void Calculate_forces (Particle *p, int number_of_particles, int cage) {
    for (int i = 0; i < number_of_particles; i++) {
        p[i].F_x = p[i].F_y = p[i].F_z = 0;
        p[i].U = 0;
    }

    for (int i = 0; i < number_of_particles; i++) {
        for (int j = 0; j < i; j++) { // Check later
            for (int a = -1; a < 2; a++) {
                for (int b = -1; b < 2; b++) {
                    for (int c = -1; c < 2; c++) {
                        double dx = p[j].x - p[i].x + a * cage;
                        double dy = p[j].y - p[i].y + b * cage;
                        double dz = p[j].z - p[i].z + c * cage;
                        double dr = sqrt(dx*dx + dy*dy + dz*dz);

                        if (dr >= cage/2.0) {
                            continue;
                        }

                        double F = - 48 * pow(dr, -13) + 24 * pow(dr, -7);
                        p[i].F_x += F * dx / dr;
                        p[i].F_y += F * dy / dr;
                        p[i].F_z += F * dz / dr;

                        p[j].F_x -= F * dx / dr;
                        p[j].F_y -= F * dy / dr;
                        p[j].F_z -= F * dz / dr;

                        // Lennard-Jones potential
                        p[i].U += 4 * (pow(dr, -12) - pow(dr, -6));
                        p[j].U += 4 * (pow(dr, -12) - pow(dr, -6));
                    }
                }
            }
        }
    }
}

void Generate_particles (Particle *p, int number_of_particles, double max_vel, double dt, int cage) {
    PRNG generator;
    initGenerator(generator);

    for (int i = 0; i < number_of_particles; i++) {
        p[i].mass = 1;

        p[i].x = p[i].x_0 = 1 * (i / (cage * cage) % (cage));
        p[i].y = p[i].y_0 = 1 * (i / cage % cage);
        p[i].z = p[i].z_0 = 1 * (i % cage);
        /*p[i].x = getRandomFloat(generator, 0, 100);
        p[i].y = getRandomFloat(generator, 0, 100);
        p[i].z = getRandomFloat(generator, 0, 100);*/

        if (i % int(cbrt(number_of_particles)) == int(cbrt(number_of_particles))-1) {
            double p_x = 0;
            double p_y = 0;
            double p_z = 0;
            for (int j = 0; j < i; j++) {
                p_x += p[j].mass * p[j].v_x;
                p_y += p[j].mass * p[j].v_y;
                p_z += p[j].mass * p[j].v_z;
            }
            p[i].v_x = - p_x / p[i].mass;
            p[i].v_y = - p_y / p[i].mass;
            p[i].v_z = - p_z / p[i].mass;
        } else {
            p[i].v_x = getRandomFloat(generator, -max_vel, max_vel);
            p[i].v_y = getRandomFloat(generator, -max_vel, max_vel);
            p[i].v_z = getRandomFloat(generator, -max_vel, max_vel);
        }

        p[i].x_prev = p[i].x - p[i].v_x * dt;
        p[i].y_prev = p[i].y - p[i].v_y * dt;
        p[i].z_prev = p[i].z - p[i].v_z * dt;

        p[i].K = (p[i].v_x*p[i].v_x + p[i].v_y*p[i].v_y + p[i].v_z*p[i].v_z) / 2;
    }
}

double Kinetic_energy (Particle *p, int number_of_particles) {
    double sum = 0;
    for (int i = 0; i < number_of_particles; i++) {
        sum += p[i].K;
    }
    return sum;
}

double Potential_energy (Particle *p, int number_of_particles) {
    double sum = 0;
    for (int i = 0; i < number_of_particles; i++) {
        sum += p[i].U;
    }
    return sum / 2;
}

double Momentum (Particle *p, int number_of_particles) {
    double p_x = 0;
    double p_y = 0;
    double p_z = 0;
    for (int i = 0; i < number_of_particles; i++) {
        p_x += p[i].mass * p[i].v_x;
        p_y += p[i].mass * p[i].v_y;
        p_z += p[i].mass * p[i].v_z;
    }
    return sqrt(p_x*p_x + p_y*p_y + p_z*p_z);
}

void Maxwell_dist (Particle *p, int *prob, int number_of_particles, double max_v, int num_segm) {
    for (int i = 0; i < number_of_particles; i++) {
        //prob[(int)floor(abs(p[i].v_x) / (max_v / num_segm))] += 1;
        int N = int(floor((p[i].v_x + max_v) / (2 * max_v / num_segm)));
        if (N < num_segm && N >=0) {
            prob[N] += 1;
        }
    }
}

double Mean_sq_r (Particle *p, int number_of_particles) {
    double mean_sq_r = 0;
    for (int i = 0; i < number_of_particles; i++) {
        mean_sq_r += (p[i].x - p[i].x_0) * (p[i].x - p[i].x_0) +
                     (p[i].y - p[i].y_0) * (p[i].y - p[i].y_0) +
                     (p[i].z - p[i].z_0) * (p[i].z - p[i].z_0);
    }
    mean_sq_r /= number_of_particles;
    return mean_sq_r;
}

int main() {
    int number_of_particles = 216; // Only cubic!
    double max_velocity = 15;
    int cage = 6;
    double dt = 0.001;
    double T = 1;
    int num_segm = 100;
    int *prob = new int[num_segm];
    for (int i = 0; i < num_segm; i++) {
        prob[i] = 0;
    }
    double mean_kin = 0;

    ofstream coords("coords.txt");
    coords << "x" << "\t" << "y" << "\t" << "z" << "\t" << "radius" << endl;
    ofstream nrg("energies.csv");
    nrg << "Kinetic energy"  << "\t" << "Potential energy" << "\t" << "Full energy" << endl;
    ofstream maxw("Distribution.csv");
    maxw << "v^2" << "\t" << "N" << "\t" << "lnN" << endl;
    ofstream r_sq("mean_sq_r.csv");
    r_sq << "t" << "\t" << "r^2" << endl;

    Particle *particles = new Particle[number_of_particles];
    Generate_particles(particles, number_of_particles, max_velocity, dt, cage);
    Calculate_forces(particles, number_of_particles, cage);

    for (double t = 0; t < T; t+=dt) {
        // Writing energies
        double K = Kinetic_energy(particles, number_of_particles);
        double P = Potential_energy(particles, number_of_particles);
        nrg << K << "\t" << P << "\t" << K + P << endl;
        mean_kin += K;
        //cout << Momentum(particles, number_of_particles) << endl;

        coords << number_of_particles << endl << endl;
        // Writing coordinates
        for (int i = 0; i < number_of_particles; i++) {
            coords << fixed << particles[i].x << ' ' << particles[i].y << ' ' << particles[i].z << ' ' << 0.2
                   << endl;
        }
        // Changing coordinates, calculating intermediate velocities
        Update_coords(particles, number_of_particles, cage, dt);

        // Calculating forces
        Calculate_forces(particles, number_of_particles, cage);

        // Calculating new velocities
        Update_velocity(particles, number_of_particles, dt);

        // Collecting data
        Maxwell_dist(particles, prob, number_of_particles,2*max_velocity, num_segm);
        r_sq << t << "\t" << Mean_sq_r(particles, number_of_particles) << endl;
    }
    // Writing data about maxwell distribution
    for (int i = 0; i < num_segm; i++) {
        maxw << pow(max_velocity * (2*i + 1)/num_segm - max_velocity, 2)
             << "\t" << prob[i] << "\t" << log(prob[i]) << endl;
    }

    // Writing mean square velocity
    mean_kin /= T/dt * number_of_particles;
    cout << mean_kin;


    coords.close();
    nrg.close();
    maxw.close();
    r_sq.close();

    delete[] prob;
    delete[] particles;
}
