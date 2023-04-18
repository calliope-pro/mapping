#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <time.h>

#include "libs/slide_min.h"

namespace fs = std::filesystem;
using libs::SlideMin;
using std::cout;
using std::endl;
using std::optional;
using std::string;
using std::unordered_map;

using Minimizer = int;
using RefId = string;
using ReadId = string;
enum class Strand : unsigned short
{
  Forward = 0,
  Reverse = 1
};
enum class Twist : unsigned short
{
  Forward = 0,
  Reverse = 1
};

const std::nullopt_t null = std::nullopt;

int main()
{
  class Minimap
  {
  private:
    const fs::path read_seq_file_path_;
    const fs::path ref_seq_file_path_;
    fs::path mapping_result_output_file_path_;
    const short window_size_;
    const short k_mer_size_;
    const int k_mer_hash_mod_;

    unordered_map<RefId, string> ref_seq_dict;
    unordered_map<ReadId, string> read_seq_dict;
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
          k_mer_hash_mod_((1 << ((k_mer_size - 1) * 2)) - 1)
    {
      time_t now;
      time(&now);
      char buf[sizeof "mapping_result_illumina_2011-10-08T07:07:09.tsv"];
      strftime(buf, sizeof buf, "mapping_result_illumina_%FT%T.tsv", gmtime(&now));
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
          tmp_ref_id = line.substr(0, idx);
        }
        else
        {
          seq += line;
        }
      }
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
      while (!ifs.eof())
      {
        string read_id;
        std::getline(ifs, read_id);
        string read_array;
        std::getline(ifs, read_array);
        std::getline(ifs, _);
        std::getline(ifs, _);
        read_seq_dict[read_id] = read_array;
      }
    }

    void run()
    {
      parse_ref_seq_file();
      parse_read_seq_file();
    }
  };

  const bool IS_TEST = false;
  const fs::path CURRENT_FOLDER_DIR = fs::current_path();

  Minimap minimap(
      CURRENT_FOLDER_DIR / "homework" / "test" / "read.fastq",
      CURRENT_FOLDER_DIR / "homework" / "test" / "ref.fasta",
      10,
      40);
  minimap.run();

  return 0;
}
