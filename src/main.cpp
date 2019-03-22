#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cstdio>
#include <ctime>
#include <typeinfo>

#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

int lucas = 0;

struct Position {
  int id;
  double x;
  double y;
  double vel;

  double s;
  double d;
};

/*
    returns a vector (len=3) of positions of cars in each lane (position 0 is cars in lane 0, etc.)
*/
vector <vector <Position>> putPositionInlane(vector<Position> closer_than_60){
  vector <vector <Position>> lanes1_to_3;

  vector <Position>lane0;
  vector <Position>lane1;
  vector <Position>lane2;

  double d;
  for(int i=0; i < closer_than_60.size(); i++){

    d = closer_than_60[i].d;

    if (d > 1.0 && d < 3.0){  
      lane0.push_back(closer_than_60[i]);
    }
    // are they in lane 1
     else if (d > 5.0 && d < 7.0){      
      lane1.push_back(closer_than_60[i]);
    }
    // are they in lane 2
    else if (d > 9.0 && d < 11.0){
      lane2.push_back(closer_than_60[i]);
    }
  }

  lanes1_to_3.push_back(lane0);
  lanes1_to_3.push_back(lane1);
  lanes1_to_3.push_back(lane2);

  return lanes1_to_3;
}

/*
    simple heuristic of which lane is safe to switch to
*/
bool isLaneSafe(vector<Position>alane, Position car_pos){

    double car_vel = car_pos.vel;
    double car_s   = car_pos.s;
    bool flag = true;

    for(int i=0; i < alane.size(); i++){
      // behind and faster
      if( (alane[i].s < car_s) && alane[i].vel > car_vel ){  
        flag = false;      
      } 
      
      //in front and slower
      if( (alane[i].s > car_s) && alane[i].vel < car_vel ){
        flag = false;         
      } 
    }

    return flag;
}

/*
  In each lane shift option, we'll define whether or not it is safe explicitly prior to running this functions
  to be on the safe side, isleftsafe and isrightsafe both start false
*/
void laneShiftProcessing(int lane, bool & isleftsafe, bool & isrightsafe, vector<Position> closer_than_60, Position car_pos) {
  // make positions lane 0,1,2 vectors
  // 
  vector <vector <Position>> allLanePositions = putPositionInlane(closer_than_60); 
  vector <Position>lane0 = allLanePositions[0];
  vector <Position>lane1 = allLanePositions[1];
  vector <Position>lane2 = allLanePositions[2];

  double car_vel = car_pos.vel;
  double car_s   = car_pos.s;

  bool flag  = false; 
  bool flag2 = false;

  if (lane == 0){
    isleftsafe = false;
    // calculate if lane1 is safe
    flag = isLaneSafe(lane1, car_pos);

    if (flag){ isrightsafe = true; }

    // see if there's a car in front within 60 going slower or behind us within x going faster  
  } else if (lane == 1){
    // calculate if lane0 and lane2 is safe
    flag  = isLaneSafe(lane0, car_pos);
    flag2 = isLaneSafe(lane2, car_pos);

    if (flag)  isleftsafe = true;
    if (flag2) isrightsafe = true;
  } else if (lane == 2){
    isrightsafe = false;
    // calculate if lane1 is safe
    flag  = isLaneSafe(lane1, car_pos);
    if (flag)  isleftsafe = true;
  }

}

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

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;

  long long start = 0; 
  int lane = 1;
  double ref_vel = 0.0;
  double spacing = .02;         
     
  bool firstlanechange = false;  

  while (getline(in_map_, line)) {
  	istringstream iss(line);
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

   
  
  h.onMessage([&firstlanechange, &start, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &ref_vel, &spacing](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
          // Under here are values           
        	// My car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

            Position car_pos;
            car_pos.id = 42;
            car_pos.x = car_x;
            car_pos.y = car_y;
            car_pos.vel = car_speed;
            car_pos.s = car_s;
            car_pos.d = car_d;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];
          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	vector <vector <double>> sensor_fusion = j[1]["sensor_fusion"];

            /* Some analysis
                    Goal here is to get a list of the vehicles within 60m of us
                    Wait, how close are they and how fast are they going for path planning we assume 0 accel on their part
             */
            vector<Position> closer_than_60;
            for( int i=0; i < sensor_fusion.size(); i++){
                
                Position p;
                p.id  = sensor_fusion[i][0];
                p.x   = sensor_fusion[i][1];
                p.y   = sensor_fusion[i][2];
                p.vel = sqrt(sensor_fusion[i][3] * sensor_fusion[i][3] + sensor_fusion[i][4] * sensor_fusion[i][4] );              
                p.d = sensor_fusion[i][6];
                double dist_to_us = distance(p.x, p.y, car_x, car_y);
                if (dist_to_us < 40){
                  closer_than_60.push_back(p);
                } 
                //cout<< "d of car id="<<p.id << " is " << p.d << endl;                 
            }

             	

            int prev_size = previous_path_x.size();

            if (prev_size > 0){
              car_s = end_path_s;
            }

            /*
              Below we check if we are too close to the car in front of us
            */
            bool too_close = false;
            start += 1;

            //https://stackoverflow.com/questions/1819189/what-range-of-values-can-integer-types-store-in-c
            // yeah probably won't get an overflow; however, we never want one
            if(start > 100000000) start -= 10000000;
            bool is_left_shift_safe  = false;
            bool is_right_shift_safe = false;

            for ( int i=0; i < sensor_fusion.size(); i++){
              float d = sensor_fusion[i][6];

              // is the car in our lane?
              if(d < 2 + 4 * lane + 2 && d > 2 + 4 * lane - 2 ){
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double check_speed = sqrt( vx*vx + vy*vy );
                double check_car_s = sensor_fusion[i][5];

                //if we're using previous points we can project the s out in time (next timestep)
                check_car_s += ((double)prev_size * .02 * check_speed);

                // is the external vehicle in front of us 
                if (check_car_s > car_s && ((check_car_s - car_s)) < 30){

                  // at this point, know we're almost in danger need to change lanes or slow down
                  
                  //the logic  (1st 3 ifs) we want to do is try to do then lane change, and if neither are safe
                  // then we'll decelerate until either it is safe to change lanes or to keep in lane

                  laneShiftProcessing( lane, is_left_shift_safe, is_right_shift_safe, closer_than_60, car_pos );

                  // Caution below
                  if (lane == 1 && (!firstlanechange || start > 70) && is_right_shift_safe){
                    
                    firstlanechange = true;

                    //decide if lane 0 or lane 2 is safer then do the safer one
                    lane = 2;                    
                    start = 0;
                  } else if (lane == 1 && (!firstlanechange || start > 35) && is_left_shift_safe){
                    //decide if safe to change into lane 1
                    firstlanechange = true;
                    lane = 0;
                    start = 0;
                    //cout << "Got here variables not modified\n";
                  }else if (lane == 2 && start > 35 && is_left_shift_safe){
                    //decide if safe to change into lane 1
                    lane = 1;
                    start = 0;
                    //cout << "Got here variables not modified\n";
                  } else if (lane == 0 && start > 35 && is_right_shift_safe){
                    //decide if safe to change into lane 1 (middle)
                    lane = 1;
                    start = 0;
                    //cout << "Got here variables not modified\n";
                  }  else{
                    // if we're going faster than them
                    if(car_speed > check_speed){
                      ref_vel -= .7;                    
                    } 
                    /*
                    else if( ref_vel < 49.5){
                      ref_vel += .224;
                    }
                    */
                  }
                  //cout << lane << endl; 
                }
              }//ends if car is in our lane

            }
            //cout<< "start  "<< start << "\t" << lane << endl;

            if(ref_vel > 49.5){
              ref_vel -= .224;
             // cout<< "Decelerationg "<< lane << endl;
            } else if( ref_vel < 49.5){
              ref_vel += .224;
            }
            

            vector <double> ptsx;
            vector <double> ptsy;

            // reference x,y yaw status
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            //double ref_vel = car_speed;            
            bool more_than_3_prev = false;

            if (prev_size < 3){
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);
              
              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);
            } else {

              more_than_3_prev = true;
              //cout<< "had 2 or more prev_points" << prev_size << endl;
              ref_x = previous_path_x[prev_size - 1];
              ref_y = previous_path_y[prev_size - 1];

              double ref_x_prev = previous_path_x[prev_size - 2];
              double ref_y_prev = previous_path_y[prev_size - 2];
              ref_yaw = atan2( (ref_y - ref_y_prev),(ref_x - ref_x_prev) );
              
              ptsx.push_back(ref_x_prev);
              ptsx.push_back(ref_x);
              
              ptsy.push_back(ref_y_prev);
              ptsy.push_back(ref_y);

            }

            /* 
            instead of the 5m spacing, they're spaced more
            .5 * .02 sec. ==> means .5m per .02sec?             
            We're on around 29:50
            */

            vector <double> next_wp0 = getXY( car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );
            vector <double> next_wp1 = getXY( car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );
            vector <double> next_wp2 = getXY( car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y );

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            // ref_x is the previous x_position, basically .02 seconds ago
            /*
                Coordinate transformation
            */
            for( int i=0;i<ptsx.size(); i++){              
              double shiftx = ptsx[i] - ref_x;
              double shifty = ptsy[i] - ref_y;

              ptsx[i] = (shiftx * cos(-1 * ref_yaw) - sin(-1 * ref_yaw) * shifty);
              ptsy[i] = (shifty * cos(-1 * ref_yaw) + sin(-1 * ref_yaw) * shiftx);

            }

            // at this point we have 2 points behind current position and 30,60,90 ahead
            tk::spline s;

            // think of ptsx and ptsy as anchor points (to smooth the curve?),
            // anchor points around which points made
            s.set_points(ptsx, ptsy);
             
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            vector<double> next_x_vals;
            vector<double> next_y_vals;

            for ( int i = 0; i < previous_path_x.size(); i++){
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }

            // code to make sure the car goes the right velocity

            double target_x = 30.0;
            double target_y = s(target_x);
            double target_dist = sqrt( (target_x * target_x) + (target_y * target_y));
            double x_add_on = 0;   

            for( int i = 1; i <= 50 - previous_path_x.size(); i++){

              double N = (target_dist / (spacing * ref_vel/2.24));
              double x_point = x_add_on + (target_x) / N ;
              double y_point = s(x_point);

              x_add_on = x_point;

              double x_ref = x_point;
              double y_ref = y_point;

              x_point = (x_ref * cos( ref_yaw) - y_ref*sin(ref_yaw));
              y_point = (x_ref * sin( ref_yaw) + y_ref*cos(ref_yaw));

              x_point += ref_x;
              y_point += ref_y;

              next_x_vals.push_back(x_point);
              next_y_vals.push_back(y_point);

            }

            //car eats up past points
            json msgJson;    

            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

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
