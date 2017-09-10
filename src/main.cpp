#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

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
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

int main() {
  uWS::Hub h;


  // MPC is initialized here!
  MPC mpc;
  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];
		  v *= 0.44704; // from MPH to m/sec

		  double steer_value_input = j[1]["steering_angle"];
		  double throttle_value_input = j[1]["throttle"];
		  steer_value_input *= deg2rad(25); // deg to rnd

          //
          // 1) We get the data in the "world" space so we need to transform them to the car space and then we work with them
          //   1.a) find the ptsx,ptsy vectors from the car point of view
          Eigen::VectorXd ptsx_car(ptsx.size());
          Eigen::VectorXd ptsy_car(ptsy.size());
		  double cos_psi = cos(psi);
		  double sin_psi = sin(psi);
		  double dx, dy;
		  for (int i=0; i<ptsx.size(); i++){
		    dx = (ptsx[i]-px);
		    dy = (ptsy[i]-py);
		    ptsx_car(i) = dx*cos_psi + dy*sin_psi;
		    ptsy_car(i) = -dx*sin_psi + dy*cos_psi;
		  }

          //
          //   1.b) find the 3ed order polynomial coeffs from the car point of view
          Eigen::VectorXd coeffs = polyfit(ptsx_car, ptsy_car, 3);

          //
          //   1.c) find the errors in the current car position (0,0)
          double cte = polyeval(coeffs, 0) - 0.0;
          double epsi = -atan(coeffs[1]);

          //
          // 2) Calculates steeering angle and throttle using MPC
          //   2.a) current state vector (x,y,psi,v,cte,epsi) from the car point of view.
          //        we will calc the car position after 100 mSec delay in order to deal with the latency of 100 mSec

		  Eigen::VectorXd state(6);
		  double latency = 0.1;
		  double Lf = 2.67;
		  double x_dly = (0.0 + v * latency);
		  double y_dly = 0.0;
		  double psi_dly = 0.0 + v * steer_value_input / Lf * latency;
		  double v_dly = 0.0 + v + throttle_value_input * latency;
		  double cte_dly = cte + (v * sin(epsi) * latency);
		  double epsi_dly = epsi + v * steer_value_input / Lf * latency;
          state << x_dly, y_dly, psi_dly, v_dly, cte_dly, epsi_dly;

		  //
		  //   2.b) use the MPC to find the next state
		  vector<double> next_state = mpc.Solve(state, coeffs);
          double steer_value = (-next_state[6])/deg2rad(25);
          double throttle_value = next_state[7];

          json msgJson;
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //
		  // 3) Display the MPC predicted trajectory and the reference line
		  //    3.a) Green line: MPC predicted trajectory
          vector<double> mpc_x_vals = mpc.result_x;
          vector<double> mpc_y_vals = mpc.result_y;
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

		  //
          //   3.b) Yellow line: Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;
		  for (int i=0; i<ptsx_car.size(); i++){
		    next_x_vals.push_back(ptsx_car(i));
		    next_y_vals.push_back(ptsy_car(i));
		  }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds((int)(latency*1000)));
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
