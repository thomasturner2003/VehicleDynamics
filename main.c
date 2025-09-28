#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/////////////////////////////////////// Parameters setup ///////////////////////////////////////
typedef struct{
    float pi;
    float inf;
}Constants;


typedef struct{
    float g;
    float mu_road;
}Environmental;


typedef struct{
    float mass;          
    float xcg;          // measured from the front axle         
    float hcg;   
    float yaw_inertia;
}Inertial;


typedef struct{
    float wheelbase;    
    float track_f;      
    float track_r;
}Geometry ;


typedef struct{
    float mu_long;
    float mu_lat;
    float sensitivity;
    float rear_stiffness;
    float front_stiffness;
}Tyres;


typedef struct{
    float roll_stiffness_f;
    float roll_stiffness_r;
    float roll_height_f;
    float roll_height_r;
}Roll;


typedef struct{
    float timestep;
}Solver;


typedef struct{
    Constants constants;
    Environmental environment;
    Inertial inertial;
    Geometry geometry;
    Tyres tyres;
    Roll roll;
    Solver solver;
}Param;

/////////////////////////////////////// State setup ///////////////////////////////////////

typedef struct{
    float x;
    float y;
    float curvature;
    float yaw;
    float time;
}Position;


typedef struct{
    float v_long;
    float v_lat;
    float yaw_rate;
    float roll_angle;
    float roll_rate;
}Dynamics;


typedef struct{
    float F_fr;
    float F_fl;
    float F_rr;
    float F_rl;
    float F_long_fr;
    float F_long_fl;
    float F_long_rr;
    float F_long_rl;
    float F_lat_fr;
    float F_lat_fl;
    float F_lat_rr;
    float F_lat_rl;
}Loads;


typedef struct{
    Position position;
    Dynamics dynamics;
    Loads loads;
}State;


/////////////////////////////////////// Functions ///////////////////////////////////////

void set_param(Param * param){
    // All units SI (radians) unless stated otherwise

    // constants
    param->constants.pi = 3.14;
    param->constants.inf = 99999;

    // environment
    param->environment.g = 9.81;
    param->environment.mu_road = 1;

    // geometry
    param->geometry.wheelbase = 1.5;
    param->geometry.track_f = 0.8 * param->geometry.wheelbase;
    param->geometry.track_r = param->geometry.track_f;

    // inertial
    param->inertial.hcg = 0.1;
    param->inertial.mass = 300;
    param->inertial.xcg = 0.5 * param->geometry.wheelbase;
    param->inertial.yaw_inertia = 100;

    // roll
    param->roll.roll_height_f = 0.1;
    param->roll.roll_height_r = 0.1;
    param->roll.roll_stiffness_f = 300;
    param->roll.roll_stiffness_r = 300;

    // solver
    param->solver.timestep = 0.01;

    // tyres
    param->tyres.mu_long = 1.4;
    param->tyres.mu_lat = 1.6;
}


void set_state(State* state, const Param* param){
    state->position.curvature = param->constants.pi/4;
    state->position.time = 0;
    state->position.x = 0;
    state->position.y = 0;
    state->position.yaw = 0;
    state->dynamics.roll_angle = 0;
    state->dynamics.roll_rate = 0;
    state->dynamics.v_lat = 0;
    state->dynamics.v_long = 1;
    state->dynamics.yaw_rate = 0;
    state->loads.F_fl = 100;
    state->loads.F_fr = 100;
    state->loads.F_rl = 100;
    state->loads.F_rr = 100;
    state->loads.F_long_fr = 0;
    state->loads.F_long_fl = 0;
    state->loads.F_long_rr = 100;
    state->loads.F_long_rl = 100;
    state->loads.F_lat_fr = 0;
    state->loads.F_lat_fl = 0;
    state->loads.F_lat_rr = 0;
    state->loads.F_lat_rl = 0;
}


void position_update(State *state, const Param *param){
    float x_dot = state->dynamics.v_long * cosf(state->position.yaw) - state->dynamics.v_lat * sinf(state->position.yaw);
    float y_dot = state->dynamics.v_long * sinf(state->position.yaw) + state->dynamics.v_lat *cosf(state->position.yaw);
    state->position.x += x_dot * param->solver.timestep;
    state->position.y += y_dot * param->solver.timestep;
    state->position.yaw += state->dynamics.yaw_rate * param->solver.timestep;
    state->position.time += param->solver.timestep;
}


float cornering_speed(const State *state, const Param *param){
    // This ignores the impact of differential by setting force of both driven loads as equal.
    // Need to fix the tyre stiffness
    float F_lat_max_fr, F_lat_max_fl, F_lat_max_rr, F_lat_max_rl;
    if (1 - powf(state->loads.F_long_fr/(param->tyres.mu_long * state->loads.F_fr), 2) < 0){
        printf("Front right exceeds limit\n");
        F_lat_max_fr = 0;
    }
    else{
        F_lat_max_fr = param->tyres.mu_lat * state->loads.F_fr * sqrtf(1 - powf(state->loads.F_long_fr/(param->tyres.mu_long * state->loads.F_fr), 2));
    }
    if (1 - powf(state->loads.F_long_fl/(param->tyres.mu_long * state->loads.F_fl), 2) < 0) {
        printf("Front left exceeds limit\n");
        F_lat_max_fl = 0;
    }
    else{
        F_lat_max_fl = param->tyres.mu_lat * state->loads.F_fl * sqrtf(1 - powf(state->loads.F_long_fl/(param->tyres.mu_long * state->loads.F_fl), 2));
    }
    if (1 - powf(state->loads.F_long_rr/(param->tyres.mu_long * state->loads.F_rr), 2) < 0){
        printf("Rear right exceeds limit\n");
        F_lat_max_rr = 0;
    }
    else{
        F_lat_max_rr = param->tyres.mu_lat * state->loads.F_rr * sqrtf(1 - powf(state->loads.F_long_rr/(param->tyres.mu_long * state->loads.F_rr), 2));
    }
    if (1 - powf(state->loads.F_long_rl/(param->tyres.mu_long * state->loads.F_rl), 2) < 0){
        printf("Rear left exceeds limit\n");
        F_lat_max_rl = 0;
    }
    else{
         F_lat_max_rl = param->tyres.mu_lat * state->loads.F_rl * sqrtf(1 - powf(state->loads.F_long_rl/(param->tyres.mu_long * state->loads.F_rl), 2));
    }
    float F_lat_capacity = F_lat_max_fr + F_lat_max_fl + F_lat_max_rr + F_lat_max_rl;
    if (state->position.curvature != 0){
        return sqrtf(F_lat_capacity/(param->inertial.mass * state->position.curvature));
    }
    else{
        return param->constants.inf;
    }
}


void wheel_loads(State *state, const Param *param, float a_long, float a_lat){
    // a_long: longitudinal acceleration, acceleration +ve, decceleration -ve
    // a_lat: lattitudal acceleration, turning right +ve, turning left -ve

    float F_static_f = (param->geometry.wheelbase - param->inertial.xcg) * param->environment.g * param->inertial.mass/param->geometry.wheelbase;
    float F_static_r = (param->inertial.xcg * param->environment.g * param->inertial.mass)/param->geometry.wheelbase;

    float F_long_f = - param->inertial.mass * a_long * param->inertial.hcg / param->geometry.wheelbase;
    float F_long_r = - F_long_f;

    float k_total = (param->roll.roll_stiffness_f + param->roll.roll_stiffness_r);
    float F_sprung_f = param->roll.roll_stiffness_f * param->inertial.mass * a_lat * param->inertial.hcg / (k_total * param->geometry.track_f);
    float F_sprung_r = param->roll.roll_stiffness_r * param->inertial.mass * a_lat * param->inertial.hcg / (k_total * param->geometry.track_r);
    float F_geom_f = param->inertial.mass * a_lat * param->roll.roll_height_f / param->geometry.track_f;
    float F_geom_r = param->inertial.mass * a_lat * param->roll.roll_height_r / param->geometry.track_r;
    float F_lat_f = F_geom_f + F_sprung_f;
    float F_lat_r = F_geom_r + F_sprung_r;

    state->loads.F_fl = (F_static_f + F_long_f - F_lat_f)/2;
    state->loads.F_fr = (F_static_f + F_long_f + F_lat_f)/2;
    state->loads.F_rl = (F_static_r + F_long_r - F_lat_r)/2;
    state->loads.F_rr = (F_static_r + F_long_r + F_lat_r)/2;
}


void drive(State *state, const Param *param){
    // selects the lateral and longitudinal forces that maximise cornering speeds
    
    // idea is to increase longitudinal force from 0 -> calculate forces -> calculate accelerations-> wheel loads -> check if you exceed the max cornering speed
}


int main() {
    State state;
    Param param;
    set_param(&param);
    set_state(&state, &param);
    state.position.yaw = 3.14/2;
    state.dynamics.v_long = 1;
    state.loads.F_long_rl = 300;
    state.loads.F_long_fl = 300;
    for(int i=0;i<100;i++){
        position_update(&state, &param);
        printf("t=%lf x=%lf y=%lf\n",state.position.time, state.position.x, state.position.y);
        cornering_speed(&state, &param);
    }
    float speed = cornering_speed(&state, &param);
    printf("%lf", speed);
    return 0;
}