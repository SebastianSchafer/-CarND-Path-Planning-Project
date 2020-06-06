#include <uWS/uWS.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include <math.h>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "helpers.h"
#include "json.hpp"
#include "spline.h"

#include "planning.cpp"

// for convenience
using nlohmann::json;
using std::string;
using std::vector;
using namespace std;


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  std::ifstream in_map_(map_file_.c_str(), std::ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    std::istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,
               &map_waypoints_dx,&map_waypoints_dy]
              (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
               uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      // this_thread::sleep_for(chrono::milliseconds(60));
      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values 
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side 
          //   of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;
// ========= insert own code here ======================================================

          car_speed = car_speed / MS_TO_MPH;
          vector<Lane> lanes;
          init_lanes(lanes);
          int remaining_steps = previous_path_x.size();
          populate_lanes(lanes, sensor_fusion, car_s, car_speed, TRANSITION_STEPS);

          // find best lane for path forward
          int target_lane = best_lane(lanes, car_lane(car_d));
          Lane &tl = lanes[target_lane];
          tk::spline spline;
          vector<double> path_s;
          vector<double> path_d;
          double angle;
          vector<double> px;
          vector<double> py;

          // path starts after transition steps to account for processing delay
          // when starting up the simulator, have no previous path and v=0
          if (remaining_steps < TRANSITION_STEPS) {
            angle = deg2rad(car_yaw);
            px.push_back(car_x);
            py.push_back(car_y);
            px.push_back(car_x + cos(angle)*1);
            py.push_back(car_y + sin(angle)*1);
          } else {
            double px1 = previous_path_x[TRANSITION_STEPS];
            double py1 = previous_path_y[TRANSITION_STEPS];
            double px0 = previous_path_x[TRANSITION_STEPS + 1];
            double py0 = previous_path_y[TRANSITION_STEPS + 1];
            px.push_back(px1);
            py.push_back(py1);
            px.push_back(px0);
            py.push_back(py0);
            car_x = previous_path_x[TRANSITION_STEPS];
            car_y = previous_path_y[TRANSITION_STEPS];
            angle = atan2(py0 - py1, px0 - px1);
            car_speed = dist(py0, px0, py1, px1) * RATE;
          }
          double t_s = target_s(car_s, min(car_speed, MAX_SPEED), tl.front_speed);
          vector<double> xy = getXY(t_s - fabs(tl.front_speed) / RATE * SPLINE_STEPS, lane_d(target_lane), 
                                    map_waypoints_s, map_waypoints_x, map_waypoints_y);
          px.push_back(xy[0]);
          py.push_back(xy[1]);
          xy = getXY(t_s, lane_d(target_lane), 
                    map_waypoints_s, map_waypoints_x, map_waypoints_y);
          px.push_back(xy[0]);
          py.push_back(xy[1]);
          double xi;

          for (int i=0; i<px.size(); i++) {
            xi = px[i];
            px[i] = (px[i] - car_x) * cos(angle) + (py[i] - car_y) * sin(angle);
            py[i] = (py[i] - car_y) * cos(angle) - (xi - car_x) * sin(angle);
          }

          spline.set_points(px, py);
          sample_spline_xy(px[0], py[0], car_speed, tl.front_speed, spline, next_x_vals, next_y_vals, 
                        map_waypoints_s, map_waypoints_x, map_waypoints_y);

          for (int i=0; i<next_y_vals.size(); i++) {
            xi = next_x_vals[i];
            next_x_vals[i] = car_x + next_x_vals[i] * cos(angle) - next_y_vals[i] * sin(angle);
            next_y_vals[i] = car_y + next_y_vals[i] * cos(angle) + xi * sin(angle);
          }

          px.clear();
          py.clear();
          if (fabs(fmod(car_d, LANE_WIDTH) - LANE_WIDTH/2) > 0.1 * LANE_WIDTH) {
            for (double d=0; d<remaining_steps; d++) {
              px.push_back(previous_path_x[d]);
              py.push_back(previous_path_y[d]);
            }
            next_x_vals = px;
            next_y_vals = py;
          } else if (remaining_steps > TRANSITION_STEPS) {
            for (double d=0; d<TRANSITION_STEPS; d++) {
              px.push_back(previous_path_x[d]);
              py.push_back(previous_path_y[d]);
            }
            px.insert(px.end(), next_x_vals.begin(), next_x_vals.end());
            py.insert(py.end(), next_y_vals.begin(), next_y_vals.end());
            next_x_vals = px;
            next_y_vals = py;
            }
            

// ============= provided code continues ==============================================
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }  // end "telemetry" if
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }  // end websocket if
  }); // end h.onMessage

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  
  h.run();
}