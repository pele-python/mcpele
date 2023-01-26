#ifndef _MCPELE_SUCCESS_CONTAINER_H
#define _MCPELE_SUCCESS_CONTAINER_H

#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcpele {

typedef unsigned long long ull;

class SuccessContainer {
 private:
  ull n_success;
  ull n_failures;

 public:
  SuccessContainer() : n_success(0), n_failures(0) {}

  inline void log_success(bool success) {
    if (success) {
      ++n_success;
    } else {
      ++n_failures;
    }
  }

  inline ull get_n_success() { return n_success; }

  inline ull get_n_failures() { return n_failures; }

  inline ull get_n_total() { return n_success + n_failures; }

  inline double get_success_rate() {
    return static_cast<double>(n_success) / (n_success + n_failures);
  }
};

class SuccessAccumulator {
 private:
  // Use std::map to store the step taken to the success of the step determined
  // later on.
  std::map<std::string, SuccessContainer> success_map;
  std::string current_step;
  bool current_success;
  bool last_success;
  bool step_being_taken;

 public:
  // Add a new string and boolean to the map.
  // SuccessAccumulator() {}

  SuccessAccumulator()
      : current_success(false), last_success(false), step_being_taken(false) {}

  void add_step_taken(std::string step) {
    current_step = step;
    step_being_taken = true;
  }

  void add_success(bool success) {
    if (!step_being_taken) {
      throw std::runtime_error("Current step not set.");
    }
    step_being_taken = false;

    if (success_map.find(current_step) == success_map.end()) {
      success_map[current_step] = SuccessContainer();
    }

    last_success = current_success;
    current_success = success;

    success_map[current_step].log_success(success);
  }

  bool get_current_success() const { return current_success; }

  bool get_last_success() const { return last_success; }

  ull get_n_successes(std::string step) {
    return success_map[step].get_n_success();
  }

  double get_success_rate(std::string step) {
    return success_map[step].get_success_rate();
  }

  ull get_n_failures(std::string step) {
    return success_map[step].get_n_failures();
  }

  std::vector<std::string> get_step_names() {
    std::vector<std::string> steps;
    for (auto it = success_map.begin(); it != success_map.end(); ++it) {
      steps.push_back(it->first);
    }
    return steps;
  }

  // For easy access to the success rates from cython
  // otherwise we would have returned std::map
  // get steps and get success rates should
  // allow us to reconstruct the map from the python
  // side
  std::vector<double> get_success_rates() {
    std::vector<double> successes;
    for (auto it = success_map.begin(); it != success_map.end(); ++it) {
      successes.push_back(it->second.get_success_rate());
    }
    return successes;
  }

  void print_success_rates() {
    std::cout << "Step Success Rate" << std::endl;
    for (auto it = success_map.begin(); it != success_map.end(); ++it) {
      std::cout << it->first << " " << it->second.get_success_rate()
                << std::endl;
    }
    std::cout << "----------------" << std::endl;
  }
};
}  // namespace mcpele

#endif  // #ifndef _MCPELE_SUCCESS_CONTAINER_H