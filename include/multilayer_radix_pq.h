//
// Created by ndrshrzg on 10/26/17.
//

//#include <stxxl/queue>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <iterator>
#include <cassert>

#ifndef MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H
#define MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H

/// Outline:
///  Bitvector of size h (limited by C in base r?) indicating whether bucket is empty
///  Buckets (std::array) of size r pointing to
///  Blocks (stxxl::queue) containing the key, val pairs
///  As pointed out in the paper the bitvector (priority queue in paper) and buckets (disk pages)
///  will fit in internal memory

//TODO methods only reading should be constexpr

namespace multilayer_radix_pq {

    namespace internal {
        // manpen: constexpr
        static constexpr std::pair<size_t, size_t> calculateBucket(uint64_t key, uint64_t last, size_t r) {
            //Left very explicit for now to debug
            if(key == last) {
                // manpen: You can just write std::make_pair(0,0) or {0,0};
                return std::pair<size_t, size_t>(0, 0);
            }
            else{
                // manpen: Declare as much as possible const
                // manpen: Wrap __builtin_clzll as it is not portable
                const size_t index_highest_significant = (std::numeric_limits<uint64_t>::digits - __builtin_clzll(key ^ last)) - 1;

                // manpen: Pow is a floating point operation; use (size_t(1) << index)
                const size_t mask = pow(2, index_highest_significant+1) -1;
                const size_t i = floor(index_highest_significant / log2(r));
                const auto shift = static_cast<size_t>(log2(r) * i);
                const size_t j = ((key & mask) >> shift);

                return std::pair<size_t, size_t>(i, j);
            }
        };


        template<typename key_type, typename block_type>
        // manpen: You pass bucket by value (involves copying). As it is an
        // container, rather use "block_type& bucket"
        static constexpr key_type scanBucketForMinimum(block_type bucket){
            key_type reorganization_minimum = std::numeric_limits<key_type>::max();
            while (!(bucket.empty())){
                if (bucket.back().first < reorganization_minimum){
                    reorganization_minimum = bucket.back().first;
                }
                bucket.pop_back();
            }
            return reorganization_minimum;
        }

        //static constexpr auto findFirstNonEmpty(/*bucket_empty_flags_*/){};

    }; // end of namespace internal

    // manpen: "C" should NOT be a template parameter, but rather a runtime
    // value. Consider to calculate upper bound on C (std::numeric_limit<KeyType>)
    // if necessary, but allow for more restricted values during runtime
    template<typename KeyType, typename ValueType, size_t RADIX_BITS, size_t C>
    class multilayer_radix_pq {
    public:
        //limiting the queue size to avoid memory overflow
        //static constexpr auto block_size = size_t(1) << 18; // deprecated for nostxxl branch
        using key_type = KeyType;
        using value_type = ValueType;
        //using block_type = stxxl::queue<std::pair<key_type, value_type>, block_size>; // deprecated for nostxxl branch
        using block_type = std::vector<std::pair<key_type, value_type>>;

    private:
        static const size_t radix = size_t(1) << RADIX_BITS;
        static const size_t no_of_buckets_ = radix-1;
        static const int no_of_arrays_ = ceil(log2(C)/RADIX_BITS);
        key_type last_minimum_;
        bool reseeding_n_flag_;

        std::array<std::array<block_type, no_of_buckets_>, no_of_arrays_> buckets_;
        std::array<std::pair<int, std::array<bool, no_of_buckets_>>, no_of_arrays_> bucket_empty_flags_;
        block_type n_bucket_;

    public:
        explicit multilayer_radix_pq() {
            bucket_empty_flags_ = {};
            last_minimum_ = std::numeric_limits<key_type>::min();
            // manpen: false
            reseeding_n_flag_ = 0;
        };

        // manpen: val should be "const value_type&"
        void push(key_type key, value_type val) {
            // TODO replace reinterpret_cast<>(key) with encoder
            // check whether monotonicity is upheld
            assert(key >= last_minimum_);

            // if key minum last minimum exceeds range [m, m+C] push into N  bucket
            if((key - last_minimum_) > C  & !reseeding_n_flag_){
                n_bucket_.push_back(std::pair<key_type, value_type> (key, val));
                std::cout << "Pushing into N bucket." << std::endl; // output
            } else {
                // calculate bucket indices (array, bucket)
                // manpen: "const auto pos"
                std::pair<uint64_t, uint64_t> pos = internal::calculateBucket(
                        reinterpret_cast<uint64_t>(key), last_minimum_, radix);

                // update bucket empty flags for calculated bucket

                // manpen:  Please write separate assignments; use "true" rather than 1
                bucket_empty_flags_[pos.first].first = bucket_empty_flags_[pos.first].second[pos.second] = 1;

                std::cout << "Pushing into B(" << pos.first << ", " << pos.second << ")." << std::endl; // output

                // push key, value pair into calculated bucket
                // manpen: at this point it is better to write VEC.emplace_back(key, val);
                buckets_[pos.first][pos.second].push_back(std::pair<key_type, value_type>(key, val));
            }
        }


        // manpen: Check semantics of PQ. In STL pop() is void and we should stick to it
        std::pair<key_type, value_type> pop() {
            assert(!empty());

            // find position of first non empty bucket
            std::pair<int, int> pos_minimum_element_ = top();
            // if regular buckets are empty
            if(pos_minimum_element_.first == -1){
                // check if N bucket holds elements
                if (n_bucket_.empty()) {
                    std::cout << "mlrpq empty" << std::endl;
                    return;
                }

                // temporary N bucket for elements outside [m,m+C] range
                block_type temp_n_bucket;
                std::cout << "begin seeding from N bucket" << std::endl;
                // set reseeding from N flag
                reseeding_n_flag_ = 1;
                // find minimum in N bucket for reorganization
                key_type reorganization_minimum = internal::scanBucketForMinimum<key_type, block_type>(n_bucket_);
                // assign new minimum
                last_minimum_ = reorganization_minimum;
                // reorganize elements from N bucket

                for(; !n_bucket_.empty(); n_bucket_.pop_back()) {
                    // retrieve arbitrary element from current bucket
                    std::pair<key_type, value_type> temp_element = n_bucket_.back();
                    // check if element is in range [m, m+C], if not it remains in N bucket
                    if (temp_element.first < (last_minimum_ + C)) {
                        // push arbitrary element in mlrpq using the calculated minimum
                        push(temp_element.first, temp_element.second);
                    }
                    // element remains in N bucket after reseeding
                    else {
                        // push element in temporary N bucket
                        temp_n_bucket.push_back(temp_element);
                    }
                }

                // set temporary N bucket as N bucket
                n_bucket_.swap(temp_n_bucket);
                std::cout << "end reseeding from N phase" << std::endl;
                // set reseeding from N flag to false
                reseeding_n_flag_ = 0;
                // call pop again after reorganization
                pop();
            }
            // if B(0, m0) empty
            else if (pos_minimum_element_.first > 0) {
                std::cout << "begin reorganization phase" << std::endl;

                // calculate minimum of current first non empty bucket
                // manpen: Type deduction should work here, so you can omit <,>
                key_type reorganization_minimum = internal::scanBucketForMinimum<key_type, block_type>(
                    buckets_[pos_minimum_element_.first][pos_minimum_element_.second]);

                // assign new minimum
                last_minimum_ = reorganization_minimum;

                // reorganize bucket into B(0, m0) using reorganization_minimum
                while(!buckets_[pos_minimum_element_.first][pos_minimum_element_.second].empty()){
                    // retrieve arbitrary element from current bucket
                    std::pair<key_type, value_type> temp_element = buckets_[pos_minimum_element_.first][pos_minimum_element_.second].back();
                    // push arbitrary element in mlrpq using the calculated minimum
                    push(temp_element.first, temp_element.second);
                    // remove element from old bucket
                    buckets_[pos_minimum_element_.first][pos_minimum_element_.second].pop_back();
                }
                // update flags for _old_ bucket
                updateBucketEmptyFlags(pos_minimum_element_);
                std::cout << "end reorganization phase" << std::endl;
                // call pop again after reorganization
                pop();
            }
            else{
                // retrieve element from bucket
                std::pair<key_type, value_type> minimum_element_ =
                        buckets_[pos_minimum_element_.first][pos_minimum_element_.second].back();
                // pop element from bucket
                buckets_[pos_minimum_element_.first][pos_minimum_element_.second].pop_back();
                // print popped element to console
                std::cout << "Popped: (" << minimum_element_.first << ", " << minimum_element_.second.x << ")."
                          << std::endl;
                // check bucket empty flags and update if necessary
                updateBucketEmptyFlags(pos_minimum_element_);
                // return minimum element
                return minimum_element_;
            }
        }

        // update bucket empty flags if vectors were emptied
        // manpen: Should this be public?
        void updateBucketEmptyFlags(const std::pair<int, int> pos_update) {
            // check if bucket is empty
            if (buckets_[pos_update.first][pos_update.second].empty()){
                // manpen: use false, rather than 0
                bucket_empty_flags_[pos_update.first].second[pos_update.second] = 0;

                // manpen: See std::count / std::count_if
                // if bucket became empty, check if whole array is empty
                // temporary variable to count empty buckets
                int temp = 0;
                for (int j = 0; j < no_of_buckets_; j++){
                    // manpen: temp += buckets_[pos_update.first][j].empty()
                    // manpen: Possibly faster

                    // if bucket is empty increase temp variable
                    if (buckets_[pos_update.first][j].empty()){ temp += 1; }
                }

                // manpen: Even better, you don't need to count, but rather use
                // std::all_of. Can be optimised better as in most cases elements
                // near the beginning will be filled!

                // if number of empty buckets equals total number of buckets assign array flag to zero
                if(temp == no_of_buckets_) { bucket_empty_flags_[pos_update.first].first = 0;}
            }
        }

        // manpen: top should return a reference to the current top ELEMENT
        // it should be either of
        //   const value_type& top() const {}
        //   const std::pair<key_type, value_type>& top() const {}
        std::pair<int,int> top() {
            //redefine as constexpr find_first_non_empty?
            for (int i = 0; i < no_of_arrays_; i++) {
                if (bucket_empty_flags_[i].first) {
                    for (int j = 0; j <= no_of_buckets_; j++){
                        if (bucket_empty_flags_[i].second[j]) {
                            return std::pair<int, int>(i,j);
                        }
                    }
                }
            }
            return std::pair<int,int>(-1,-1);
        }

        // manpen: empty() should be REALLY fast, O(1) time
        // manpen: make it this const, i.e. bool empty() const {}
        bool empty() {
            return (top() == std::pair<int,int>(-1, -1));
        }

        //this method only exists for a unit test which is weird
        bool bucketEmpty(int k){
            return !bucket_empty_flags_[k].first;
        }

    };
}


#endif //MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H
