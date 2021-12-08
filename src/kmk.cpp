#include <algorithm>
#include <cassert>
#include <cstdint>
#include <list>
#include <memory>
#include <vector>
#include <Rcpp.h>

using Rcpp::DataFrame;
using Rcpp::Named;

// the kk heuristic, from: https://stackoverflow.com/questions/68421928/how-to-reconstruct-partitions-in-the-karmarkar-karp-heuristic-multi-way-partitio
// to do: transform to full kk?


namespace {



class Subset {
 public:
  Subset() : numbers_{}, sum_{} {}
  explicit Subset(const int number) : numbers_{number}, sum_{number} {}

  Subset(Subset&&) = default;
  Subset& operator=(Subset&&) = default;

  const std::list<int>& numbers() const { return numbers_; }

  int sum() const { return sum_; }

  void Merge(Subset other) {
    numbers_.splice(numbers_.end(), other.numbers_);
    sum_ += other.sum_;
  }

 private:
  Subset(const Subset&) = delete;
  Subset& operator=(const Subset&) = delete;

  std::list<int> numbers_;
  int sum_;
};


struct SubsetSumGreater {
  bool operator()(const Subset& a, const Subset& b) const {
    return a.sum() > b.sum();
  }
};


class Partition {
 public:
  Partition(const int number, const std::size_t k) : subsets_(k) {
    assert(k > 0);
    subsets_.front().Merge(Subset{number});
  }

  Partition(Partition&&) = default;
  Partition& operator=(Partition&&) = default;

  const std::vector<Subset>& subsets() const { return subsets_; }

  int difference() const {
    return subsets_.front().sum() - subsets_.back().sum();
  }

  void Merge(Partition other) {
    assert(subsets_.size() == other.subsets_.size());
    auto it{subsets_.begin()};
    auto other_it{other.subsets_.rbegin()};
    while (it != subsets_.end() && other_it != other.subsets_.rend()) {
      it->Merge(std::move(*other_it));
      ++it;
      ++other_it;
    }
    std::sort(subsets_.begin(), subsets_.end(), SubsetSumGreater{});
  }

 private:
  Partition(const Partition&) = delete;
  Partition& operator=(const Partition&) = delete;

  std::vector<Subset> subsets_;
};


struct DifferenceLess {
  bool operator()(const Partition& a, const Partition& b) const {
    return a.difference() < b.difference();
  }
};


Partition KarmarkarKarp(const std::vector<int>& numbers,
                        const std::size_t k) {
  assert(!numbers.empty());
  assert(k > 0);
  std::vector<Partition> heap;
  heap.reserve(numbers.size());
  for (const int number : numbers) {
    heap.emplace_back(number, k);
  }
  std::make_heap(heap.begin(), heap.end(), DifferenceLess{});
  while (heap.size() > 1) {
    std::pop_heap(heap.begin(), heap.end(), DifferenceLess{});
    Partition partition{std::move(heap.back())};
    heap.pop_back();
    std::pop_heap(heap.begin(), heap.end(), DifferenceLess{});
    heap.back().Merge(std::move(partition));
    std::push_heap(heap.begin(), heap.end(), DifferenceLess{});
  }
  return std::move(heap.front());
}

}  // namespace



// [[Rcpp::export]]
DataFrame kmk_partition(const std::vector<int>& nrs, const int k)
{
	int j=0, g=1;
	std::vector<int> ogr(nrs.size());
	std::vector<int> onr(nrs.size());
	
	const Partition kk = KarmarkarKarp(nrs, k);	
	
	for(auto ss{kk.subsets().begin()}; ss != kk.subsets().end(); ++ss)
	{
		for(auto it{(*ss).numbers().begin()}; it != (*ss).numbers().end(); ++it)
		{
			onr[j] = *it;
			ogr[j++] = g;
		}
		g++;		
	}

	return DataFrame::create(Named("size")=onr, Named("group")=ogr);
}



