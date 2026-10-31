// Included only by temporary sources built for --profile.
#pragma once

#include <array>
#include <chrono>
#include <ostream>
#include <vector>

namespace nep_benchmark
{
enum Stage { setup, neighbors, descriptor, radial, angular, bec, electrostatics, zbl, stage_count };

struct Batch {
  std::array<double, stage_count> seconds = {};
  int calls = 0;
};

struct State {
  bool enabled = false;
  Batch current;
  std::vector<Batch> batches;
};

inline State& state()
{
  static State value;
  return value;
}

class Call
{
public:
  Call() : enabled_(state().enabled)
  {
    if (enabled_) {
      last_ = Clock::now();
    }
  }

  ~Call()
  {
    if (enabled_) {
      start(setup);
      ++state().current.calls;
    }
  }

  void start(Stage next)
  {
    if (enabled_) {
      const auto now = Clock::now();
      state().current.seconds[stage_] += std::chrono::duration<double>(now - last_).count();
      last_ = now;
      stage_ = next;
    }
  }

private:
  using Clock = std::chrono::steady_clock;
  bool enabled_;
  Stage stage_ = setup;
  Clock::time_point last_;
};

inline void begin_batch()
{
  state().current = Batch();
  state().enabled = true;
}

inline void end_batch()
{
  state().enabled = false;
  state().batches.push_back(state().current);
}

inline void write(std::ostream& out)
{
  out << ",\"stage_seconds\":[";
  for (std::size_t i = 0; i < state().batches.size(); ++i) {
    if (i) {
      out << ',';
    }
    out << '[';
    for (int stage = 0; stage < stage_count; ++stage) {
      if (stage) {
        out << ',';
      }
      out << state().batches[i].seconds[stage];
    }
    out << ']';
  }
  out << "],\"stage_calls\":[";
  for (std::size_t i = 0; i < state().batches.size(); ++i) {
    if (i) {
      out << ',';
    }
    out << state().batches[i].calls;
  }
  out << ']';
}
} // namespace nep_benchmark
