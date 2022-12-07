#include <cstddef>
#include <map>
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
    } 
    else {
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

public:
  // Add a new string and boolean to the map.
  SuccessAccumulator() {}

  void add_step_taken(std::string step) { current_step = step; }

  void add_success(bool success) {

    if (success_map.find(current_step) == success_map.end()) {
      success_map[current_step] = SuccessContainer();
    }

    last_success = current_success;
    current_success = success;

    success_map[current_step].log_success(success);
  }

  bool get_current_success() { return current_success; }

  bool get_last_success() { return last_success; }

  ull get_n_success(std::string step) {
    return success_map[step].get_n_success();
  }

  double get_success_rate(std::string step) {
    return success_map[step].get_success_rate();
  }

  ull get_n_failures(std::string step) {
    return success_map[step].get_n_failures();
  }

  std::vector<std::string> get_steps() {
    std::vector<std::string> steps;
    for (auto it = success_map.begin(); it != success_map.end(); ++it) {
      steps.push_back(it->first);
    }
    return steps;
  }

};
} // namespace mcpele