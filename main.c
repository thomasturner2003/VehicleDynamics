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
    float xcg;                   
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


float cornering_speed(State *state, const Param *param){
    // This ignores the impact of differential by setting force of both driven loads as equal.
    // Need to fix the tyre stiffness
    float F_lat_max_fr = param->tyres.mu_lat * state->loads.F_fr * sqrtf(1 - powf(state->loads.F_long_fr/(param->tyres.mu_long * state->loads.F_fr), 2));
    float F_lat_max_fl = param->tyres.mu_lat * state->loads.F_fl * sqrtf(1 - powf(state->loads.F_long_fl/(param->tyres.mu_long * state->loads.F_fl), 2));
    float F_lat_max_rr = param->tyres.mu_lat * state->loads.F_rr * sqrtf(1 - powf(state->loads.F_long_rr/(param->tyres.mu_long * state->loads.F_rr), 2));
    float F_lat_max_rl = param->tyres.mu_lat * state->loads.F_rl * sqrtf(1 - powf(state->loads.F_long_rl/(param->tyres.mu_long * state->loads.F_rl), 2));
    float F_lat_capacity = F_lat_max_fr + F_lat_max_fl + F_lat_max_rr + F_lat_max_rl;
    if (state->position.curvature != 0){
        return sqrtf(F_lat_capacity/(param->inertial.mass * state->position.curvature));
    }
    else{
        return param->constants.inf;
    }
}

void load_transfer(State *state, const Param *param, float a_long, float a_lat, bool counter_clockwise){
    // longitudinal
    float delta_F_long_fr = param->inertial.mass * a_long * param->inertial.hcg/2;
    float delta_F_long_fl = delta_F_long_fr;
    float delta_F_long_rr = - delta_F_long_fr;
    float delta_F_long_rl = - delta_F_long_fr;

    // lattitudal
    float k_total = (param->roll.roll_stiffness_f + param->roll.roll_stiffness_r);
    float delta_F_sprung_f = param->roll.roll_stiffness_f * param->inertial.mass * a_lat * param->inertial.hcg / (k_total * param->geometry.track_f);
    float delta_F_sprung_r = param->roll.roll_stiffness_r * param->inertial.mass * a_lat * param->inertial.hcg / (k_total * param->geometry.track_r);
    float delta_F_geom_f = param->inertial.mass * a_lat * param->roll.roll_height_f;
    float delta_F_geom_r = param->inertial.mass * a_lat * param->roll.roll_height_r;
    float delta_F_f = delta_F_geom_f + delta_F_sprung_f;
    float delta_F_r = delta_F_geom_r + delta_F_sprung_r;

    // adjust the state loads based on h
}

int main() {
    State state;
    Param param;
    set_param(&param);
    set_state(&state, &param);
    state.position.yaw = 3.14/2;
    state.dynamics.v_long = 1;
    for(int i=0;i<10;i++){
        position_update(&state, &param);
        printf("t=%lf x=%lf y=%lf\n",state.position.time, state.position.x, state.position.y);
    }
    float speed = cornering_speed(&state, &param);
    printf("%lf", speed);
    return 0;
}