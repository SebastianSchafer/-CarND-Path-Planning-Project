

#include <iostream>
#include "helpers.h"
#include "spline.h"
#include <math.h>


using namespace std;

using std::vector;

const int LANE_WIDTH = 4;
const int NUMBER_OF_LANES = 3;
const double TRACK_LENGTH = 6945.554;
const double MS_TO_MPH = 2.23694;
const double MAX_SPEED = 49.5 / MS_TO_MPH;
const double MAX_ACC = 10;
const double MAX_JERK = 10;
const double RATE = 50;
const double SENSOR_RANGE = 150;
const double FRONT_CLEARANCE = 15; // time/speed based for real app
const double REAR_CLEARANCE = 15;
const double PATH_LENGTH = 2.5; // length of new path in s
const int SPLINE_STEPS = 5;
const int TRANSITION_STEPS = 10; 

struct Lane {
	bool blocked;
	double front_clearance;
	double front_speed;
	double rear_clearance;
	double rear_speed;
	double cost;
};


double dist(double x1, double y1, double x2, double y2) {
  return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

//returns int of lane, 0-based for indexing
int car_lane(double ego_d) {
	return (int)(ego_d / LANE_WIDTH);
}

double lane_d(int lane) {
	return 2 + lane*LANE_WIDTH;
}

void clear_lane(Lane &lane) {
	lane.blocked = false;
	lane.front_clearance = SENSOR_RANGE;
	lane.front_speed = MAX_SPEED;
	lane.rear_clearance = SENSOR_RANGE;
	lane.rear_speed = MAX_SPEED;
	lane.cost = 100;
}

// reset lanes before populating with data
void init_lanes(vector<Lane> &lanes) {
	for (int i=0; i < NUMBER_OF_LANES; i++) {
		Lane lane;
		clear_lane(lane);
		lanes.push_back(lane);
	}
}

// populate lanes with sensor data; first clear lanes, then find closest vehicles
void populate_lanes(vector<Lane> &lanes, vector<vector<double>> sensor_fusion, double ego_s, double ego_v, int steps_ahead) {
	for (auto lane : lanes) {
		clear_lane(lane);
	}
	for(auto car : sensor_fusion) {
		if (car[6] < 0) continue; // simulator initially has garbage neg data
		int l_id = car_lane(car[6]);
		double car_v = sqrt(pow(car[3],2) + pow(car[4],2)) / MS_TO_MPH; 
		double car_sr = car[5] - ego_s + (car_v) * steps_ahead / RATE;

		if (fabs(car_sr) <= max(FRONT_CLEARANCE, REAR_CLEARANCE)) lanes[l_id].blocked = true;
		if (car_sr < SENSOR_RANGE && car_sr >= 0) {
			lanes[l_id].front_clearance = min(lanes[l_id].front_clearance, fabs(car_sr));
			lanes[l_id].front_speed = min(lanes[l_id].front_speed,  MAX_SPEED - 0.5*(MAX_SPEED-car_v)*FRONT_CLEARANCE/car_sr);
		}
		if (car_sr <  FRONT_CLEARANCE && car_sr >= 0) {
			lanes[l_id].front_speed = min(lanes[l_id].front_speed, car_v);
		}
		// if car behind
		if (fabs(car_sr) < SENSOR_RANGE && car_sr <= 0) {
			lanes[l_id].rear_clearance = min(lanes[l_id].rear_clearance, car_sr);
			lanes[l_id].rear_speed = max(lanes[l_id].rear_speed, car_v);
		}
		// For warpping around track length; c++ does not even have a decent mod operator:
		if (fabs(ego_s - TRACK_LENGTH) <= SENSOR_RANGE) {
			double car_sr = fmod(fmod(car[5] - ego_s, TRACK_LENGTH) + TRACK_LENGTH, TRACK_LENGTH) + (car_v) * steps_ahead / RATE;
			lanes[l_id].front_clearance = min(lanes[l_id].front_clearance, fabs(car_sr));
			lanes[l_id].blocked = true;
		}
	}
}


// Various cost functions:

// lane speed cost
double cost_v(double v) {
	return 1 * (1 - min(v / MAX_SPEED, MAX_SPEED));
}

// lane change cost; favor center lane, followed by left lane
double cost_lc(int target_lane, int ego_lane) {
	vector<double> lane_cost = {1e-2, -1e-2, 2e-2};
	if (target_lane == ego_lane) return 0;
	return lane_cost[target_lane];

}
// lane occupancy cost
double cost_occ(Lane lane, int i, int ego_lane) {
	if (i == ego_lane) {
		return cost_v(lane.front_speed);
	} 
	if (lane.blocked == true) {
		return 1e2 + cost_v(lane.front_speed);
	}
	return cost_v(lane.front_speed);
}

// find lane with lowest cost; return index of target lane
int best_lane(vector<Lane> &lanes, int ego_lane) {
	vector<double> costs;
	costs.clear();
	for (int i = 0; i < NUMBER_OF_LANES; i++) {
		if (fabs(i - ego_lane) > 1) {
			lanes[i].cost = 1e3;
			costs.push_back(lanes[i].cost);
			continue;
		}
		lanes[i].cost = cost_occ(lanes[i], i, ego_lane) + cost_lc(i, ego_lane);
		costs.push_back(lanes[i].cost);
	}
	for ( auto i : costs) cout << i << " ";
	cout << endl;
	int minel = std::distance(costs.begin(), std::min_element(costs.begin(), costs.end()));
	return minel;
}

// Find trget distance to travel; used for appr. spline length and for seeting speed on spline sampling
double target_s(double ego_s, double ego_v, double target_v, double t_total = PATH_LENGTH) {
	ego_v = min(ego_v, MAX_SPEED);
	if (fabs(ego_v / target_v - 1) < 0.05) {
		return ego_s + target_v * t_total;
	}
	double af = 1;
	double acc_t = min(fabs(target_v - ego_v)/MAX_ACC, t_total); // shortest time to ge to target velocity; not always best choice
	if (1/RATE == t_total) af = 1;
	if (target_v > ego_v){
		return ego_s + (t_total - acc_t) * target_v + acc_t * ego_v + af*0.5 * pow(acc_t,2) * MAX_ACC;
	} else{
		return ego_s + (t_total - acc_t) * target_v + acc_t * ego_v - 0.5 * pow(acc_t,2) * MAX_ACC; // ToDo: dec. slower to maximize speed?
	}
}

// sample spline for max v and a; assume s only for simplicity
void sample_spline(double start_s, double ego_v, double target_v, tk::spline spline,
					vector<double> &next_x_vals, 
					vector<double> &next_y_vals, 
					const vector<double> &maps_s, 
					const vector<double> &maps_x, 
					const vector<double> &maps_y) {
	vector<double> xy;
	double next_s = start_s;
	xy = getXY(start_s,spline(start_s), maps_s, maps_x, maps_y);
	next_x_vals.push_back(xy[0]);
	next_y_vals.push_back(xy[1]);
	for (int dt = 1; dt/RATE <= PATH_LENGTH; dt++) {
		start_s = next_s;
		next_s = target_s(start_s, ego_v, target_v, 1/RATE);
		ego_v = (next_s - start_s) * RATE; // update ego_v for next iter
		xy = getXY(next_s,spline(next_s), maps_s, maps_x, maps_y);
		next_x_vals.push_back(xy[0]);
		next_y_vals.push_back(xy[1]);
	}

}

// sample spline in xy coords
void sample_spline_xy(double start_x, double start_y, double ego_v, double target_v, tk::spline spline,
					vector<double> &next_x_vals, 
					vector<double> &next_y_vals, 
					const vector<double> &maps_s, 
					const vector<double> &maps_x, 
					const vector<double> &maps_y) {
	double next_dist;
	double x0 = start_x;
	double y0 = spline(start_x);
	double next_x = x0;
	double next_y = y0;
	double x_step;
	double curr_angle;
	double prev_angle = 0;
	double prev_v;
	double accel;
	next_x_vals.push_back(x0);
	next_y_vals.push_back(y0);
	for (int dt = 1; dt/RATE <= PATH_LENGTH; dt++) {
		start_x = next_x;
		start_y = next_y;
		next_dist = target_s(0, ego_v, target_v, 1/RATE);
		x_step = next_dist;
		next_x =  start_x + x_step;
		next_y = spline(next_x);
		if (fabs((next_y-start_y)/x_step) > 0.15) {
			next_dist = target_s(0, ego_v, ego_v, 1/RATE);
		}
		while (dist(start_x, start_y, next_x, next_y) >= next_dist) {
			x_step *= 0.999;
			next_x =  start_x + x_step;
			next_y = spline(next_x);
			}
		ego_v = dist(start_x, start_y, next_x, next_y) * RATE; // update ego_v for next iter
		next_x_vals.push_back(next_x + x0);
		next_y_vals.push_back(next_y + y0);
	}

}
