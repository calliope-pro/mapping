#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>
#include <time.h>

#include "libs/slide_min.h"

namespace fs = std::filesystem;
using libs::SlideMin;
using sc = std::chrono::system_clock;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::tuple;
using std::unordered_map;
using std::vector;

using Minimizer = unsigned long long;
using RefId = string;
using ReadId = string;
enum class Strand : unsigned short
{
  Forward = 0,
  Reverse = 1
};

const unordered_map<char, short> NUCLEOTIDE_TO_INTEGER_MAPPING({{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}});
const std::nullopt_t null = std::nullopt;

class Minimap
{
private:
  const fs::path read_seq_file_path_;
  const fs::path ref_seq_file_path_;
  fs::path mapping_result_output_file_path_;
  const short window_size_;
  const short k_mer_size_;
  const unsigned long long k_mer_hash_mod_;

  unordered_map<RefId, string> ref_seq_dict;
  unordered_map<ReadId, string> read_seq_dict;
  unordered_map<Minimizer, vector<tuple<RefId, libs::Index, Strand>>> ref_minimizer_dict;
  // int processes_;

public:
  Minimap(
      const fs::path read_seq_file_path,
      const fs::path ref_seq_file_path,
      const short window_size,
      const short k_mer_size,
      optional<fs::path> mapping_result_output_file_path = null)
      : read_seq_file_path_(read_seq_file_path),
        ref_seq_file_path_(ref_seq_file_path),
        window_size_(window_size),
        k_mer_size_(k_mer_size),
        k_mer_hash_mod_((1ULL << ((k_mer_size - 1) * 2)) - 1)
  {
    time_t now;
    time(&now);
    char buf[sizeof "mapping_result_PacBio_2011-10-08T07:07:09.tsv"];
    strftime(buf, sizeof buf, "mapping_result_PacBio_%FT%T.tsv", gmtime(&now));
    mapping_result_output_file_path_ = (mapping_result_output_file_path.has_value()) ? mapping_result_output_file_path.value() : buf;
  }

  void parse_ref_seq_file()
  {
    std::ifstream ifs(ref_seq_file_path_);
    if (ifs.fail())
    {
      std::cerr << "Failed to open file." << std::endl;
      std::exit(1);
    }
    string tmp_ref_id;
    string seq = "";
    while (!ifs.eof())
    {
      string line;
      std::getline(ifs, line);
      if (line[0] == '>')
      {
        if (seq != "")
        {
          ref_seq_dict[tmp_ref_id] = seq;
          seq = "";
        }
        const short idx = line.find_first_of(' ');
        tmp_ref_id = line.substr(1, idx-1);
      }
      else
      {
        seq += line;
      }
    }
    ifs.close();
    ref_seq_dict[tmp_ref_id] = seq;
  }

  void parse_read_seq_file()
  {
    std::ifstream ifs(read_seq_file_path_);
    if (ifs.fail())
    {
      std::cerr << "Failed to open file." << std::endl;
      std::exit(1);
    }
    string tmp_ref_id;
    string seq = "";
    string _;
    while (true)
    {
      string read_id;
      std::getline(ifs, read_id);
      string read_array;
      std::getline(ifs, read_array);
      std::getline(ifs, _);
      std::getline(ifs, _);
      if (ifs.eof())
      {
        ifs.close();
        break;
      }
      read_seq_dict[read_id.substr(1)] = read_array;
    }
  }

  void sketch_minimizers(string &seq, std::queue<tuple<libs::Value, libs::Index, Strand>> &queue)
  {
    std::set<libs::Index> minimizers_set;
    const int seq_length = seq.size();
    const string k_mer_init = seq.substr(0, k_mer_size_);
    SlideMin k_mer_hashes_strand_0(seq2hash(k_mer_init), window_size_);
    SlideMin k_mer_hashes_strand_1(seq2hash(k_mer_init, Strand::Reverse), window_size_);

    for (short i = 0; i < window_size_ - 1; i++)
    {
      const Minimizer &last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.last_value();
      const Minimizer &last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.last_value();
      Minimizer new_k_mer_hash_strand_0_value;
      Minimizer new_k_mer_hash_strand_1_value;
      calc_hash(
          last_k_mer_hash_strand_0_value,
          last_k_mer_hash_strand_1_value,
          seq[k_mer_size_ + i],
          new_k_mer_hash_strand_0_value,
          new_k_mer_hash_strand_1_value);
      k_mer_hashes_strand_0.add(new_k_mer_hash_strand_0_value);
      k_mer_hashes_strand_1.add(new_k_mer_hash_strand_1_value);
    }
    const auto &[k_mer_hash_strand_0_index,
          k_mer_hash_strand_0_value] = k_mer_hashes_strand_0.min();
    const auto &[k_mer_hash_strand_1_index,
          k_mer_hash_strand_1_value] = k_mer_hashes_strand_1.min();
    if (k_mer_hash_strand_0_value < k_mer_hash_strand_1_value)
    {
      minimizers_set.insert(k_mer_hash_strand_0_index);
      queue.push(tuple(k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, Strand::Forward));
    }
    if (k_mer_hash_strand_0_value > k_mer_hash_strand_1_value)
    {
      minimizers_set.insert(k_mer_hash_strand_1_index);
      queue.push(tuple(k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, Strand::Reverse));
    }

    for (int j = 0; j < seq_length - k_mer_size_ - window_size_ + 1; j++)
    {
      const Minimizer &last_k_mer_hash_strand_0_value = k_mer_hashes_strand_0.last_value();
      const Minimizer &last_k_mer_hash_strand_1_value = k_mer_hashes_strand_1.last_value();
      Minimizer new_k_mer_hash_strand_0_value;
      Minimizer new_k_mer_hash_strand_1_value;
      calc_hash(
          last_k_mer_hash_strand_0_value,
          last_k_mer_hash_strand_1_value,
          seq[k_mer_size_ + window_size_ + j - 1],
          new_k_mer_hash_strand_0_value,
          new_k_mer_hash_strand_1_value);
      const auto &[k_mer_hash_strand_0_index,
            k_mer_hash_strand_0_value] = k_mer_hashes_strand_0.add(new_k_mer_hash_strand_0_value);
      const auto &[k_mer_hash_strand_1_index,
            k_mer_hash_strand_1_value] = k_mer_hashes_strand_1.add(new_k_mer_hash_strand_1_value);
      if (
          k_mer_hash_strand_0_value < k_mer_hash_strand_1_value && minimizers_set.count(k_mer_hash_strand_0_index) == 0)
      {
        minimizers_set.insert(k_mer_hash_strand_0_index);
        queue.push(tuple(k_mer_hash_strand_0_value, k_mer_hash_strand_0_index, Strand::Forward));
      }
      if (
          k_mer_hash_strand_0_value > k_mer_hash_strand_1_value && minimizers_set.count(k_mer_hash_strand_1_index) == 0)
      {
        minimizers_set.insert(k_mer_hash_strand_1_index);
        queue.push(tuple(k_mer_hash_strand_1_value, k_mer_hash_strand_1_index, Strand::Reverse));
      }
    }
  }

  string generate_analyzed_read_seq_output(ReadId read_id, string read_seq)
  {
    unsigned int read_seq_length = read_seq.size();
    unordered_map<RefId, unordered_map<libs::Index, unsigned short>> hits;
    std::queue<tuple<libs::Value, libs::Index, Strand>> read_minimizer_queue;
    sketch_minimizers(read_seq, read_minimizer_queue);
    while (!read_minimizer_queue.empty())
    {
      auto &[read_minimizer, read_pos, read_strand] = read_minimizer_queue.front();
      for (auto &[ref_id, ref_pos, ref_strand] : ref_minimizer_dict[read_minimizer])
      {
        if (read_strand == ref_strand)
        {
          hits[ref_id][ref_pos - read_pos + 1]++;
        }
        else
        {
          hits[ref_id][-(ref_pos - (read_seq_length - k_mer_size_ - read_pos) + 1)]++;
        }
      }
      read_minimizer_queue.pop();
    }

    short max_cnt;
    string output;
    for (auto &[ref_id, counter] : hits)
    {
      const auto &most_common_pair = std::max_element(
          counter.begin(), counter.end(), [](const auto &x, const auto &y)
          { return (x.second < y.second); });
      if (most_common_pair->second > max_cnt)
      {
        max_cnt = most_common_pair->second;
        output = read_id + '\t' + ref_id + '\t' + std::to_string(std::abs(most_common_pair->first)) + '\t' + "+-"[most_common_pair->first < 0];
      }
    }
    return output;
  }

  void calc_hash(
      Minimizer last_k_mer_hash_strand_0_value,
      Minimizer last_k_mer_hash_strand_1_value,
      char next_nucleotide,
      Minimizer &k_mer_hash_strand_0_value,
      Minimizer &k_mer_hash_strand_1_value)
  {
    k_mer_hash_strand_0_value = ((last_k_mer_hash_strand_0_value & k_mer_hash_mod_) << 2) | seq2hash(string(1, next_nucleotide), Strand::Forward);
    k_mer_hash_strand_1_value = (last_k_mer_hash_strand_1_value >> 2) | (seq2hash(string(1, next_nucleotide), Strand::Reverse) << (2 * (k_mer_size_ - 1)));
  }

  unsigned long long seq2hash(string seq, Strand strand = Strand::Forward)
  {
    unsigned long long res = 0;
    if (strand == Strand::Forward)
    {
      for (char &base : seq)
      {
        res = (res << 2) | NUCLEOTIDE_TO_INTEGER_MAPPING.at(base);
      }
    }
    else
    {
      for (short idx = seq.size(); idx-- > 0;)
      {
        char &base = seq[idx];
        res = (res << 2) | (3 - NUCLEOTIDE_TO_INTEGER_MAPPING.at(base));
      }
    }
    return res;
  }

  void run()
  {
    parse_ref_seq_file();
    parse_read_seq_file();

    std::queue<tuple<libs::Value, libs::Index, Strand>> ref_minimizer_queue;
    for (auto &ref_seq_pair : ref_seq_dict)
    {
      sketch_minimizers(ref_seq_pair.second, ref_minimizer_queue);
      while (!ref_minimizer_queue.empty())
      {
        auto &[minimizer, ref_pos, strand] = ref_minimizer_queue.front();
        ref_minimizer_dict[minimizer].push_back(tuple(ref_seq_pair.first, ref_pos, strand));
        ref_minimizer_queue.pop();
      }
    }

    std::ofstream ofs(mapping_result_output_file_path_);
    for (auto &[read_id, read_seq] : read_seq_dict)
    {
      string output = generate_analyzed_read_seq_output(read_id, read_seq);
      ofs << output << '\n';
    }
    ofs.close();
  }
};

int main()
{
  const bool IS_TEST = false;
  const fs::path CURRENT_FOLDER_DIR = fs::current_path();
  if (IS_TEST)
  {
    Minimap minimap(
        CURRENT_FOLDER_DIR / "homework" / "test" / "read.fastq",
        CURRENT_FOLDER_DIR / "homework" / "test" / "ref.fasta",
        10,
        15);
    minimap.run();
  }
  else
  {
    Minimap minimap(
        CURRENT_FOLDER_DIR / "SE11" / "PacBio_SE11.fastq",
        CURRENT_FOLDER_DIR / "SE11" / "ref_SE11.fasta",
        10,
        15);
    minimap.run();
  }

  return 0;
}
