//
// Created by Andreas Herzog on 12/07/17.
// andreas.herzog@stud.uni-frankfurt.de
//

#include <stxxl/queue>
//#include <queue>
#include <array>
#include <cmath>
#include <limits>
#include <iterator>
#include <cassert>
#include <climits>

#ifndef MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H
#define MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H

#define LIMITMEMORY
//#define COMPILERAGNOSTIC

namespace internal {

    template<typename KeyType>
    class calculations{
    public:
        using key_type = KeyType;

        static constexpr size_t calculateIndexHighestSignificantBit(key_type key){
            size_t index = 0;
            for (size_t i = sizeof(key) * CHAR_BIT; i--; )
            {
                if ((key >> i) == 1){
                    return index;
                }
                index++;
            }
            return index;
        }
    };
}



namespace multilayer_radix_pq {

    template<typename KeyType, typename ValueType, size_t RADIX_BITS>
    class multilayer_radix_pq {

    public:
        using key_type = KeyType;
        using value_type = ValueType;
#ifdef LIMITMEMORY
        static constexpr auto block_size = size_t(1) << 18;
        using block_type = stxxl::queue<std::pair<key_type, value_type>, block_size>;
#else
        using block_type = stxxl::queue<std::pair<key_type, value_type>>;
#endif


    private:
        static const size_t radix_ = size_t(1) << RADIX_BITS; // manpen: constexpr? (also for later instances)
        static const size_t no_of_buckets_ = radix_;

        // estimate C as upper limit of key_type
        //static const size_t C = std::numeric_limits<key_type>::max();
        static const size_t C_ = size_t(1) << 20;
        static const size_t no_of_arrays_ = ceil(log2(C_)/RADIX_BITS) + 1;
        key_type last_minimum_;
        key_type N_bucket_minimum_;
        std::pair<int64_t, int64_t> current_minimum_index_;
        std::array<std::array<block_type, no_of_buckets_>, no_of_arrays_> buckets_;
        std::array<std::array<uint64_t, no_of_buckets_>, no_of_arrays_> bucket_minimum_;
        std::array<std::pair<bool, std::array<bool, no_of_buckets_>>, no_of_arrays_> bucket_empty_flags_;
        block_type n_bucket_;
        bool reseeding_n_flag_;
        bool first_push_flag_;


    public:
        // manpen: Use initialisation list
        // explicit multilayer_radix_pq() :
        //   bucket_empty_flags_({}),
        //   last_minimum_(std::numeric_limits<key_type>::min()),
        //   ...
        // {}
        explicit multilayer_radix_pq() {
            bucket_empty_flags_ = {};
            last_minimum_ = std::numeric_limits<key_type>::min();
            reseeding_n_flag_ = false;
            first_push_flag_ = true;
            initializeBucketMinima(no_of_arrays_, no_of_buckets_);
            N_bucket_minimum_ = std::numeric_limits<key_type>::max();
            current_minimum_index_.first = -1;
            current_minimum_index_.second = -1;

        };


        // manpen: What about copy/move constructors and assignments?

        ~multilayer_radix_pq() = default;


    private:
        // manpen: write size_t /* buck */ to avoid unused warning or remove param
        void initializeBucketMinima(size_t arr, size_t buck){
            // initialize array with bucket minima to be maximum
            for (int i=0; i < arr; i++){
                bucket_minimum_[i].fill(std::numeric_limits<key_type>::max());
            }
        }

        void reseedFromNBucket() {
            // re initialize bucket minima
            initializeBucketMinima(no_of_arrays_, no_of_buckets_);
            // temporary N bucket for elements outside [m,m+C] range
            block_type temp_n_bucket;
            // set reseeding from N flag
            reseeding_n_flag_ = 1;
            // retrieve minimum in N bucket for reorganization
            //key_type reorganization_minimum = bucket_minimum_[pos_minimum_element_.first][pos_minimum_element_.second];
            const key_type reorganization_minimum = N_bucket_minimum_;
            // assign new minimum
            last_minimum_ = reorganization_minimum;
            // reset N_bucket_minimum_
            N_bucket_minimum_ = std::numeric_limits<key_type>::max();
            // reorganize elements from N bucket
            for (; !n_bucket_.empty(); n_bucket_.pop()) {
                // retrieve arbitrary element from current bucket
                std::pair<key_type, value_type> temp_element = n_bucket_.front();
                // check if element is in range [m, m+C], if not it remains in N bucket

                if (temp_element.first < (last_minimum_ + C_)) {
                    // push arbitrary element in mlrpq using the calculated minimum
                    push(temp_element.first, temp_element.second);

                } else {
                    // manpen: If we have the invariant that no element is larger than
                    // current-min + C, and the queue structure can accommodate
                    // at least an interval of  size C, we do NOT need to push
                    // elements in temp_n_bucket. While it's a nice feature, I
                    // would strongly suggest to make a template parameter which
                    // by default forbids this behaviour and then use an assertion to check
                    // that no index is too large to overflow N (otherwise the whole IO analysis
                    // falls apart). Include this parameter also in this if
                    // statement to allow the compiler to optimise it away in
                    // case overflowing is disallowed


                    // element remains in N bucket after reseeding
                    // push element in temporary N bucket
                    temp_n_bucket.push(temp_element);
                    if (temp_element.first < N_bucket_minimum_) {
                        N_bucket_minimum_ = temp_element.first;
                    };
                }
            }

            // set temporary N bucket as N bucket
            n_bucket_.swap(temp_n_bucket);
            // set reseeding from N flag to false
            reseeding_n_flag_ = 0;
            // call pop again after reorganization
        }

        // manpen: It somehow seems unnatural to use a pair rather than
        // two arguments here. Won't impact performance though

        // manpen: You really want to implement it using word-wise parallelism
        void updateBucketEmptyFlags(const std::pair<int64_t, int64_t> pos_update) {
            // check if bucket is empty
            if (buckets_[pos_update.first][pos_update.second].empty()){
                bucket_empty_flags_[pos_update.first].second[pos_update.second] = 0;

                // manpen: As already discussed you do not want to count, but
                // rather find a non-empty bucket; so simply break the loop
                // on the first non-empty bucket

                // if bucket became empty, check if whole array is empty
                // temporary variable to count empty buckets
                int temp = 0;
                for (int j = 0; j < no_of_buckets_; j++){
                    // if bucket is empty increase temp variable
                    if (buckets_[pos_update.first][j].empty()){ temp += 1; }
                }
                // if number of empty buckets equals total number of buckets assign array flag to zero
                if(temp == no_of_buckets_) { bucket_empty_flags_[pos_update.first].first = 0;}
            }
        }


        const std::pair <size_t, size_t> calculateBucket(uint64_t key, uint64_t last) const {
            if(key == last) {
                return {0, 0};
            }
            else{
                // manpen: It better to encapsulate this switching logic into
                // its own function, e.g. by defining "calculateIndexHighestSignificantBit"
                // either as hand-written or by using an intrinsic. Also observe
                // that STXXL gives a lot of preprocessor macros to detect which
                // compiler you are using

#ifdef COMPILERAGNOSTIC
                const size_t index_highest_significant =
                        (std::numeric_limits<uint64_t>::digits - internal::calculations<key_type>::calculateIndexHighestSignificantBit(key ^ last)) - 1;
#else
                const size_t index_highest_significant = (std::numeric_limits<uint64_t>::digits - __builtin_clzll(key ^ last)) - 1;
#endif

                // manpen: log2(radix_) = RADIX_WIDTH -- but as float!
                // Please remove any float-point arith; dividing by integers is
                // automatically floor!
                const size_t mask = (size_t(1) << (index_highest_significant+1))-1;
                const size_t i = floor(index_highest_significant / log2(radix_));
                const auto shift = static_cast<size_t>(log2(radix_) * i);
                const size_t j = ((key & mask) >> shift);
                return {i, j};
            }
        }


        // manpen: const this-pointer!
        std::pair<int64_t, int64_t> calculateMinimumIndex() {
            // manpen: You really want to implement that more efficient using
            // word-wise bit parallelism (e.g. keep a 32 or 64 empty flags
            // in a word and find high sig bit and stuff ;))
            for (int i = 0; i < no_of_arrays_; i++) {
                if (bucket_empty_flags_[i].first) {
                    for (int j = 0; j <= no_of_buckets_; j++){
                        if (bucket_empty_flags_[i].second[j]) {
                            return std::pair<int64_t, int64_t>(i,j);
                        }
                    }
                }
            }
            return {-1,-1};
        }


        void updateCurrentMinimumIndex(){
            // manpen: why not directly: current_minimum_index_ = calculateMinimumIndex()?

            std::pair<int64_t, int64_t> calculated_index = calculateMinimumIndex();
            current_minimum_index_.first = calculated_index.first;
            current_minimum_index_.second = calculated_index.second;
        }


    public:
        void push(key_type key, const value_type& val) {
            // TODO replace reinterpret_cast<>(key) with encoder
            // manpen: Yep, please to so!

            // check whether monotonicity is upheld
            assert(key >= 0);
            assert(key >= last_minimum_);

            // if key minimum last minimum exceeds range [m, m+C] push into N  bucket
            if((key - last_minimum_) > C_  & !reseeding_n_flag_){
                // manpen: does {key, val} not work here?
                n_bucket_.push(std::pair<key_type, value_type> (key, val));
                if(key < N_bucket_minimum_) {N_bucket_minimum_ = key;};
                first_push_flag_ = false;

                // manpen: Add a return here and pull up the else-scope; removes
                // one level of indention
            }
            else{
                // calculate bucket indices (array, bucket); also reinterpret_cast
                // is the wrong one. If any, use static_cast; see e.g.,
                // https://stackoverflow.com/questions/332030/when-should-static-cast-dynamic-cast-const-cast-and-reinterpret-cast-be-used
                const auto pos = calculateBucket(
                        reinterpret_cast<uint64_t>(key), last_minimum_);

                // update bucket empty flags for calculated bucket
                bucket_empty_flags_[pos.first].first = true;
                bucket_empty_flags_[pos.first].second[pos.second] = true;

                // manpen: You are allowed to use more than one line for an if ;)
                // update bucket minimum if necessary
                if (key < bucket_minimum_[pos.first][pos.second]) {bucket_minimum_[pos.first][pos.second] = key;};

                // push key, value pair into calculated bucket
                buckets_[pos.first][pos.second].push({key, val});

                // manpen: You should avoid ifs that need to be checked for every
                // element but are only relevant for the first/last item. In this
                // case you can simply initialise current_min_index with num_lim::max,
                // or am I missing something? Also try to get rid of the reseeding case.

                // set the current minimum index to calculated bucket of first ever push
                if (first_push_flag_) {
                    current_minimum_index_.first = pos.first;
                    current_minimum_index_.second = pos.second;
                }
                // when reseeding from N, current minimum index will always be (0, 0)
                else if (reseeding_n_flag_){
                    current_minimum_index_.first = 0;
                    current_minimum_index_.second = 0;
                }
                // update current minimum index without calling the expensive update function
                else{
                    if (pos.first <= current_minimum_index_.first){
                        if (pos.second < current_minimum_index_.second){
                            current_minimum_index_.first = pos.first;
                            current_minimum_index_.second = pos.second;
                        }
                    }
                }
                first_push_flag_ = false;
            }
        }


        void pop() {
            assert(!empty());
            // find position of first non empty bucket
            const auto pos_minimum_element_ = current_minimum_index_;
            // if regular buckets are empty
            if(pos_minimum_element_.first == -1){
                reseedFromNBucket();
                pop();
                // manpen: return and remove else
            }

            // if B(0, m0) empty
            else if (pos_minimum_element_.first > 0) {
                // get minimum of current first non empty bucket
                const key_type reorganization_minimum = bucket_minimum_[pos_minimum_element_.first][pos_minimum_element_.second];

                // manpen: define a shortcut for the bucket you're interested in
                // auto& min_bucket = buckets_[pos_minimum_element_.first][pos_minimum_element_.second];

                // assign new minimum
                last_minimum_ = reorganization_minimum;

                // re-initialize bucket minima up to current bucket
                initializeBucketMinima(pos_minimum_element_.first, no_of_buckets_);

                // reorganize bucket into B(0, m0) using reorganization_minimum
                for(    ;
                        !buckets_[pos_minimum_element_.first][pos_minimum_element_.second].empty();
                        buckets_[pos_minimum_element_.first][pos_minimum_element_.second].pop()) {

                    // manpen: const auto& temp_element = ...

                    // retrieve arbitrary element from current bucket
                    std::pair<key_type, value_type> temp_element = buckets_[pos_minimum_element_.first][pos_minimum_element_.second].front();
                    // push arbitrary element in mlrpq using the calculated minimum
                    push(temp_element.first, temp_element.second);
                }

                // update flags for _old_ bucket
                updateBucketEmptyFlags(pos_minimum_element_);

                // call pop again after reorganization
                pop();

                // manpen: return and pull out else
            }
            else{
                // retrieve element from bucket
                std::pair<key_type, value_type> minimum_element_ =
                        buckets_[pos_minimum_element_.first][pos_minimum_element_.second].front();
                // pop element from bucket
                buckets_[pos_minimum_element_.first][pos_minimum_element_.second].pop();
                // check bucket empty flags and update if necessary
                updateBucketEmptyFlags(pos_minimum_element_);
                // update current minimum index
                updateCurrentMinimumIndex();

                // manpen: here, you can remove the return for a change ;)
                return;
            }
        }


        const std::pair<key_type, value_type>& top() {
            assert(!empty());

            // manpen: && rather than & !
            // manpen: Also, it seems that !empty() suffices as it includes the first condition
            if ((current_minimum_index_ == std::pair<int64_t, int64_t> (-1, -1)) & !empty()){
                reseedFromNBucket();
                top();
            }

            //std::pair<int64_t, int64_t> pos_minimum_emelent = current_minimum_index_;
            std::pair<key_type, value_type>& minimum_element = buckets_[current_minimum_index_.first][current_minimum_index_.second].front();
            return minimum_element;
        };


        bool empty() const {
            // manpen: && rather than & !
            return (current_minimum_index_ == std::pair<int64_t,int64_t>(-1, -1)) & n_bucket_.empty();
        }

        // manpen: Please provide a possibility to query (at least) the key of the
        // smallest element without changing current_minimum_index (which happens
        // in top() in case you need to reseed). Also, this function should be pretty fast

    }; // end class multilayer_radix_pq
}


#endif //MULTILAYER_RADIX_PRIORITY_QUEUE_MULTILAYER_RADIX_PQ_H
