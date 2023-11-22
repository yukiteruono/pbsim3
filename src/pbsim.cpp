#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>

#define SUCCEEDED 1
#define FAILED 0
#define TRUE 1
#define FALSE 0
#define BUF_SIZE 10240
#define REF_ID_LEN_MAX 128
#define TRANS_ID_LEN_MAX 128
#define REF_SEQ_NUM_MAX 9999
#define REF_SEQ_LEN_MAX 1000000000
#define REF_SEQ_LEN_MIN 100
#define FASTQ_NUM_MAX 100000000
#define FASTQ_LEN_MAX 1000000
#define TEMPLATE_NUM_MAX 100000000
#define TEMPLATE_LEN_MAX 1000000
#define TRANS_LEN_MAX 1000000
#define EXISTS_LINE_FEED 1
#define DEPTH_MAX 1000
#define RATIO_MAX 1000
#define STRATEGY_WGS 1
#define STRATEGY_TRANS 2
#define STRATEGY_TEMPL 3
#define METHOD_QS 1
#define METHOD_ERR 2
#define METHOD_SAM 3
#define METHOD_SAM_REUSE 4
#define METHOD_SAM_STORE 5
#define ACCURACY_MAX 100
#define STATE_MAX 50
#define TR_RANK_MAX 1000

/////////////////////////////////////////
// Definitions of structures           //
/////////////////////////////////////////

// simulation parameters
struct sim_t {
  int set_flg[30];
  int strategy;
  int method;
  unsigned int seed;
  double depth;
  double accuracy_mean, accuracy_max, accuracy_min;
  long len_min, len_max; 
  double len_mean, len_sd; 
  long long len_quota;
  long sub_ratio, ins_ratio, del_ratio;
  double sub_rate, ins_rate, del_rate;
  long res_num;
  long long res_len_total; 
  double res_depth;
  double res_accuracy_mean, res_accuracy_sd;
  long res_len_min, res_len_max; 
  double res_len_mean, res_len_sd; 
  long res_sub_num, res_ins_num, res_del_num;
  double res_sub_rate, res_ins_rate, res_del_rate;
  char *prefix, *outfile_ref, *outfile_fq, *outfile_maf, *outfile_sam;
  char *profile_id, *profile_fq, *profile_stats;
  char *id_prefix;
  int pass_num, res_pass_num;
  double hp_del_bias;
};

// Sample reads
struct sample_t {
  char *file;
  long num;
  long long len_total;
  long len_min, len_max;
  long num_filtered;
  long long len_total_filtered;
  long len_min_filtered, len_max_filtered;
  double len_mean_filtered, len_sd_filtered;
  double accuracy_mean_filtered, accuracy_sd_filtered;
};

// Genome sequences
struct genome_t {
  char *file;
  char *seq;
  short *hp;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num_seq;
  long num;
  long hpfreq[11];
  double hp_del_bias[11];
};

// Transcript sequences
struct transcript_t {
  char *file;
  char *seq;
  short *hp;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num_seq;
  long plus_exp, minus_exp, total_exp;
  int rank_max;
  long hpfreq[11];
  double hp_del_bias[11];
};

// Template sequences
struct templ_t {
  char *file;
  char *seq;
  short *hp;
  char id[REF_ID_LEN_MAX + 1];
  long len;
  long num;
  long long len_total;
  long hpfreq[11];
  double hp_del_bias[11];
};

// Mutation
struct mut_t {
  long sub_thre[94], ins_thre[94], del_thre[94], err_thre[94];
  char *sub_nt_a, *sub_nt_t, *sub_nt_g, *sub_nt_c, *sub_nt_n, *ins_nt;
  char *err;
  char *qc, *new_qc, *tmp_qc;
  char *seq, *read_seq, *maf_seq, *maf_ref_seq;
  short *hp;
  long len;
  char seq_strand;
  long seq_left, seq_right;
  long err_num, del_num;
  long offset;
  int acc;
};

// Quality code
struct qc_t {
  char character;
  double prob;
};

// Error
struct err_t {
  char character;
};

// Quality code of HMM
struct qshmm_t {
  char *file;
  double ip[ACCURACY_MAX+1][STATE_MAX+1];
  double ep[ACCURACY_MAX+1][STATE_MAX+1][94];
  double tp[ACCURACY_MAX+1][STATE_MAX+1][STATE_MAX+1];
  int exist_hmm[ACCURACY_MAX+1];
};

// Error of HMM
struct errhmm_t {
  char *file;
  double ip[ACCURACY_MAX+1][STATE_MAX+1];
  double ep[ACCURACY_MAX+1][STATE_MAX+1][4];
  double tp[ACCURACY_MAX+1][STATE_MAX+1][STATE_MAX+1];
  int exist_hmm[ACCURACY_MAX+1];
  int state_max[ACCURACY_MAX+1];
  int acc_min;
  int acc_max;
};

/////////////////////////////////////////
// Global variances                    //
/////////////////////////////////////////

FILE *fp_filtered, *fp_stats, *fp_fq, *fp_sam, *fp_maf;
struct sim_t sim;
struct sample_t sample;
struct genome_t genome;
struct transcript_t transcript;
struct templ_t templ;
struct mut_t mut;
struct qc_t qc[94];
struct err_t err[4];
struct qshmm_t qshmm;
struct errhmm_t errhmm;
long freq_len[FASTQ_LEN_MAX + 1];
long freq_accuracy[100000 + 1];
double uni_ep[ACCURACY_MAX+1][94];

/////////////////////////////////////////
// Prototypes of functions             //
/////////////////////////////////////////

int trim(char *line);
void init_sim_res();
int set_sim_param();
int set_qshmm();
int set_errhmm();
int set_mut();
int get_sample_inf();
int get_genome_inf();
int get_genome_seq();
int get_transcript_inf();
int get_templ_inf();
int simulate_by_sample();
int simulate_by_qshmm();
int simulate_by_qshmm_trans();
int simulate_by_qshmm_templ();
int simulate_by_errhmm();
int simulate_by_errhmm_trans();
int simulate_by_errhmm_templ();
int count_digit(long num);
void print_sim_param();
void print_transcript_stats();
void print_sample_stats();
void print_templ_stats();
void print_simulation_stats();
void print_help();
void revcomp(char* str);
void revshort(short* str, long len);
long get_time_cpu();
long get_time();

/////////////////////////////////////////
// Main                                //
/////////////////////////////////////////

int main (int argc, char** argv) {
  char *tp, *tmp_buf;
  long len;
  long num;
  long ratio;
  long i, j, k;
  long rst1, rst2;
  long t1, t2;
  double prob, rate;
  long sum1, sum2;

  rst1 = get_time_cpu();
  t1 = get_time();

  memset(sim.set_flg, 0, sizeof(sim.set_flg));

  sim.seed = (unsigned int)time(NULL);

  // Variables for Option
  int opt, option_index;
  struct option long_options[] = {
    {"strategy", 1, NULL, 0},
    {"method", 1, NULL, 0},
    {"genome", 1, NULL, 0},
    {"transcript", 1, NULL, 0},
    {"prefix", 1, NULL, 0},
    {"id-prefix", 1, NULL, 0},
    {"depth", 1, NULL, 0},
    {"length-min", 1, NULL, 0},
    {"length-max", 1, NULL, 0},
    {"difference-ratio", 1, NULL, 0},
    {"seed", 1, NULL, 0},
    {"sample", 1, NULL, 0},
    {"sample-profile-id", 1, NULL, 0},
    {"accuracy-min", 1, NULL, 0},
    {"accuracy-max", 1, NULL, 0},
    {"qshmm", 1, NULL, 0},
    {"errhmm", 1, NULL, 0},
    {"length-mean", 1, NULL, 0},
    {"length-sd", 1, NULL, 0},
    {"accuracy-mean", 1, NULL, 0},
    {"pass-num", 1, NULL, 0},
    {"template", 1, NULL, 0},
    {"hp-del-bias", 1, NULL, 0},
    {0, 0, 0, 0}
  };

  // Option parsing
  option_index = 0;
  while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
    switch (opt) {
    case 0:
      sim.set_flg[option_index] = 1;

      switch (option_index) {
      case 0:
        if (strncmp(optarg,"wgs",3) == 0) {
          sim.strategy = STRATEGY_WGS;
        } else if (strncmp(optarg,"trans",5) == 0) {
          sim.strategy = STRATEGY_TRANS;
        } else if (strncmp(optarg,"templ",5) == 0) {
          sim.strategy = STRATEGY_TEMPL;
        } else {
          fprintf(stderr, "ERROR (strategy: %s): Acceptable value: wgs, trans, templ.\n", optarg);
          exit(-1);
        }
        break;

      case 1:
        if (strncmp(optarg,"qshmm",5) == 0) {
          sim.method = METHOD_QS;
        } else if (strncmp(optarg,"errhmm",6) == 0) {
          sim.method = METHOD_ERR;
        } else if (strncmp(optarg,"sample",6) == 0) {
          sim.method = METHOD_SAM;
        } else {
          fprintf(stderr, "ERROR (method: %s): Acceptable value: qshmm, errhmm, sample.\n", optarg);
          exit(-1);
        }
        break;

      case 2:
        if ((genome.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(genome.file, optarg);
        break;

      case 3:
        if ((transcript.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(transcript.file, optarg);
        break;

      case 4:
        if ((sim.prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.prefix, optarg);
        break;

      case 5:
        if ((sim.id_prefix = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.id_prefix, optarg);
        break;

      case 6:
        sim.depth = atof(optarg);
        if (sim.depth <= 0.0) {
          fprintf(stderr, "ERROR (depth: %s): Acceptable range is more than 0.\n", optarg);
          exit(-1);
        }
        break;

      case 7:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_min = atoi(optarg);
        if ((sim.len_min < 1) || (sim.len_min > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-min: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 8:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        sim.len_max = atoi(optarg);
        if ((sim.len_max < 1) || (sim.len_max > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-max: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 9:
        if ((tmp_buf = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(tmp_buf, optarg);
        num = 0;
        tp = strtok(tmp_buf, ":");
        while (num < 3) {
          if (tp == NULL) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Format is sub:ins:del.\n", optarg);
            exit(-1);
          }
          if (strlen(tp) >= 5) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          ratio = atoi(tp);
          if ((ratio < 0) || (ratio > RATIO_MAX)) {
            fprintf(stderr, "ERROR (difference-ratio: %s): Acceptable range is 0-%d.\n", optarg, RATIO_MAX);
            exit(-1);
          }
          if (num == 0) {
            sim.sub_ratio = ratio;
          } else if (num == 1) {
            sim.ins_ratio = ratio;
          } else if (num == 2) {
            sim.del_ratio = ratio;
          }
          num ++;
          tp = strtok(NULL, ":");
        }
        break;

      case 10:
        sim.seed = (unsigned int)atoi(optarg);
        break;

      case 11:
        if ((sample.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sample.file, optarg);
        break;

      case 12:
        if ((sim.profile_id = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(sim.profile_id, optarg);
        break;

      case 13:
        sim.accuracy_min = atof(optarg);
        if ((sim.accuracy_min < 0.0) || (sim.accuracy_min > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-min: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 14:
        sim.accuracy_max = atof(optarg);
        if ((sim.accuracy_max < 0.0) || (sim.accuracy_max > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-max: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 15:
        if ((qshmm.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(qshmm.file, optarg);
        break;

      case 16:
        if ((errhmm.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(errhmm.file, optarg);
        break;

      case 17:
        sim.len_mean = atof(optarg);
        if ((sim.len_mean < 1) || (sim.len_mean > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-mean: %s): Acceptable range is 1-%d.\n", optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 18:
        sim.len_sd = atof(optarg);
        if ((sim.len_sd < 0) || (sim.len_sd > FASTQ_LEN_MAX)) {
          fprintf(stderr, "ERROR (length-sd: %s): Acceptable range is 0-%d.\n",
            optarg, FASTQ_LEN_MAX);
          exit(-1);
        }
        break;

      case 19:
        sim.accuracy_mean = atof(optarg);
        if ((sim.accuracy_mean < 0.0) || (sim.accuracy_mean > 1.0)) {
          fprintf(stderr, "ERROR (accuracy-mean: %s): Acceptable range is 0.0-1.0.\n", optarg);
          exit(-1);
        }
        break;

      case 20:
        sim.pass_num = atoi(optarg);
        if (sim.pass_num < 1) {
          fprintf(stderr, "ERROR (pass_num: %s): Acceptable range is more than 1.\n", optarg);
          exit(-1);
        }
        break;

      case 21:
        if ((templ.file = (char *)malloc(strlen(optarg) + 1)) == 0) {
          fprintf(stderr, "ERROR: Cannot allocate memory.\n");
          exit(-1);
        }
        strcpy(templ.file, optarg);
        break;

      case 22:
        if (strlen(optarg) >= 8) {
          fprintf(stderr, "ERROR (hp-del-bias: %s): Acceptable range is 1-10.\n", optarg);
          exit(-1);
        }
        sim.hp_del_bias = atof(optarg);
        if ((sim.hp_del_bias < 1) || (sim.hp_del_bias > 10)) {
          fprintf(stderr, "ERROR (hp-del-bias: %s): Acceptable range is 1-10.\n", optarg);
          exit(-1);
        }
        break;

      default:
        break;
      }
    break;

    default:
      exit(-1);
      break;
    }
  }

  if (argc == 1) {
    print_help();
    exit(-1);
  }

  // Setting of simulation parameters     
  if (set_sim_param() == FAILED) {
    exit(-1);
  }
  print_sim_param();

  srand((unsigned int)sim.seed);

  // Quality code to error probability
  for (i=0; i<=93; i++) {
    qc[i].prob = pow(10, (double)i / -10);
    qc[i].character = (char)(i+33);
  }

  // Error to error character
  err[0].character = '0';
  err[1].character = '1';
  err[2].character = '2';
  err[3].character = '3';

  // Uniform error probability
  for (i=0; i<=ACCURACY_MAX; i++) {
    for (j=0; j<=93; j++) {
      uni_ep[i][j] = 0;
    }
    if (i == ACCURACY_MAX) {
      uni_ep[i][93] = 1.0;
      continue;
    }
    prob = 1.0 - i / 100.0;
    for (j=0; j<=93; j++) {
      if (prob == qc[j].prob) {
        uni_ep[i][j] = 1.0;
        break;
      } else if (prob > qc[j].prob) {
        rate = (prob - qc[j].prob) / (qc[j-1].prob - qc[j].prob);
        uni_ep[i][j-1] = rate;
        uni_ep[i][j] = 1 - rate;
        break;
      }
    }
  }

  // Sample reads
  if (sim.method == METHOD_SAM) {
    if ((fp_filtered = tmpfile()) == NULL) {
      fprintf(stderr, "ERROR: Cannot open temporary file\n");
      return FAILED;
    }
    if (get_sample_inf() == FAILED) {
      exit(-1);
    }
    print_sample_stats();
  } else if (sim.method == METHOD_SAM_STORE) {
    if ((fp_filtered = fopen(sim.profile_fq, "w+")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "w+")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if (get_sample_inf() == FAILED) {
      exit(-1);
    }
    print_sample_stats();
  } else if (sim.method == METHOD_SAM_REUSE) {
    if ((fp_filtered = fopen(sim.profile_fq, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if ((fp_stats = fopen(sim.profile_stats, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open sample_profile\n");
      return FAILED;
    }
    if (get_sample_inf() == FAILED) {
      exit(-1);
    }
    print_sample_stats();
  }

  // Quality HMM
  if (sim.method == METHOD_QS) {
    if (set_qshmm() == FAILED) {
      exit(-1);
    }
  }
  if(0){
  for(i=0; i<=ACCURACY_MAX; i++) {
    for(j=1; j<=STATE_MAX; j++) {
      for(k=0; k<=93; k++) {
        fprintf(stderr,"EP %ld : %ld : %ld : %f\n",i,j,k,qshmm.ep[i][j][k]);
      }
    }
    for(j=1; j<=STATE_MAX; j++) {
      for(k=1; k<=STATE_MAX; k++) {
        fprintf(stderr,"TP %ld : %ld : %ld : %f\n",i,j,k,qshmm.tp[i][j][k]);
      }
    }
  }
  }

  // Error HMM
  if (sim.method == METHOD_ERR) {
    if (set_errhmm() == FAILED) {
      exit(-1);
    }
  }
  if(0){
  for(i=errhmm.acc_min; i<=errhmm.acc_max; i++) {
    for(j=1; j<=errhmm.state_max[i]; j++) {
      for(k=0; k<=3; k++) {
        fprintf(stderr,"EP %ld : %ld : %ld : %f\n",i,j,k,errhmm.ep[i][j][k]);
      }
      fprintf(stderr,"\n");
    }
    for(j=1; j<=errhmm.state_max[i]; j++) {
      for(k=1; k<=errhmm.state_max[i]; k++) {
        fprintf(stderr,"TP %ld : %ld : %ld : %f\n",i,j,k,errhmm.tp[i][j][k]);
      }
    }
  }
  }

  // Set mutation parameters and varianeces
  if (set_mut() == FAILED) {
    exit(-1);
  }

  // WGS
  if (sim.strategy == STRATEGY_WGS) {

    if (get_genome_inf() == FAILED) {
      exit(-1);
    }

    if (sim.hp_del_bias == 1) {
      for(i=1; i<=10; i++) {
        genome.hp_del_bias[i] = 1;
      }
    } else {
      for(i=0; i<=10; i++) {
        genome.hpfreq[i] = 0;
      }
      for (genome.num=1; genome.num<=genome.num_seq; genome.num++) {
        if (get_genome_seq() == FAILED) {
          exit(-1);
        }
      }
      sum1 = 0;
      sum2 = 0;
      for(i=1; i<=10; i++) {
        genome.hp_del_bias[i] = 1 + (sim.hp_del_bias - 1) / 9 * (i - 1);
        sum1 += genome.hpfreq[i] * genome.hp_del_bias[i];
        sum2 += genome.hpfreq[i];
      }
      rate = (double)sum2 / sum1;
      for(i=1; i<=10; i++) {
        genome.hp_del_bias[i] *= rate;
      }
    }

    for (genome.num=1; genome.num<=genome.num_seq; genome.num++) {
      if (get_genome_seq() == FAILED) {
        exit(-1);
      }

      init_sim_res();
      sim.len_quota = (long long)(sim.depth * genome.len);

      if (sim.pass_num == 1) {
        sprintf(sim.outfile_fq, "%s_%04ld.fastq", sim.prefix, genome.num);
        if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
          fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
          return FAILED;
        }
      } else {
        sprintf(sim.outfile_sam, "%s_%04ld.sam", sim.prefix, genome.num);
        if ((fp_sam = fopen(sim.outfile_sam, "w")) == NULL) {
          fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_sam);
          return FAILED;
        }
        fprintf(fp_sam, "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n");
        fprintf(fp_sam, "@RG\tID:ffffffff\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000\tPU:%s%ld\tPM:SEQUELII\n", sim.id_prefix, genome.num);
      }

      sprintf(sim.outfile_maf, "%s_%04ld.maf", sim.prefix, genome.num);
      if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
        return FAILED;
      }

      if (sim.method == METHOD_QS) {
        if (simulate_by_qshmm() == FAILED) {
          exit(-1);
        }
      } else if (sim.method == METHOD_ERR) {
        if (simulate_by_errhmm() == FAILED) {
          exit(-1);
        }
      } else {
        if (simulate_by_sample() == FAILED) {
          exit(-1);
        }
      }

      print_simulation_stats();

      if (sim.pass_num == 1) {
        fclose(fp_fq);
      } else {
        fclose(fp_sam);
      }
      fclose(fp_maf);
    }

    if ((sim.method == METHOD_SAM_STORE) || (sim.method == METHOD_SAM_REUSE)) {
      fclose(fp_filtered);
      fclose(fp_stats);
    }
  // TRANSCRIPT
  } else if (sim.strategy == STRATEGY_TRANS) {

    if (get_transcript_inf() == FAILED) {
      exit(-1);
    }
    print_transcript_stats();

    init_sim_res();

    if (sim.pass_num == 1) {
      sprintf(sim.outfile_fq, "%s.fastq", sim.prefix);
      if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
        return FAILED;
      }
    } else {
      sprintf(sim.outfile_sam, "%s.sam", sim.prefix);
      if ((fp_sam = fopen(sim.outfile_sam, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_sam);
        return FAILED;
      }
      fprintf(fp_sam, "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n");
      fprintf(fp_sam, "@RG\tID:ffffffff\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000\tPU:%s\tPM:SEQUELII\n", sim.id_prefix);
    }

    sprintf(sim.outfile_maf, "%s.maf", sim.prefix);
    if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
      return FAILED;
    }

    if (sim.method == METHOD_QS) {
      if (simulate_by_qshmm_trans() == FAILED) {
        exit(-1);
      }
    } else if (sim.method == METHOD_ERR) {
      if (simulate_by_errhmm_trans() == FAILED) {
        exit(-1);
      }
    }

    print_simulation_stats();

    if (sim.pass_num == 1) {
      fclose(fp_fq);
    } else {
      fclose(fp_sam);
    }
    fclose(fp_maf);
  // TEMPLATE
  } else {

    if (get_templ_inf() == FAILED) {
      exit(-1);
    }
    print_templ_stats();

    init_sim_res();

    if (sim.pass_num == 1) {
      sprintf(sim.outfile_fq, "%s.fastq", sim.prefix);
      if ((fp_fq = fopen(sim.outfile_fq, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_fq);
        return FAILED;
      }
    } else {
      sprintf(sim.outfile_sam, "%s.sam", sim.prefix);
      if ((fp_sam = fopen(sim.outfile_sam, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_sam);
        return FAILED;
      }
      fprintf(fp_sam, "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n");
      fprintf(fp_sam, "@RG\tID:ffffffff\tPL:PACBIO\tDS:READTYPE=SUBREAD;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-789-500;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000\tPU:%s\tPM:SEQUELII\n", sim.id_prefix);
    }

    sprintf(sim.outfile_maf, "%s.maf", sim.prefix);
    if ((fp_maf = fopen(sim.outfile_maf, "w")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_maf);
      return FAILED;
    }

    if (sim.method == METHOD_QS) {
      if (simulate_by_qshmm_templ() == FAILED) {
        exit(-1);
      }
    } else if (sim.method == METHOD_ERR) {
      if (simulate_by_errhmm_templ() == FAILED) {
        exit(-1);
      }
    }

    print_simulation_stats();

    if (sim.pass_num == 1) {
      fclose(fp_fq);
    } else {
      fclose(fp_sam);
    }
    fclose(fp_maf);
  }

  rst2 = get_time_cpu();
  t2 = get_time();

  fprintf(stderr, ":::: System utilization ::::\n\n");
  fprintf(stderr, "CPU time(s) : %ld\n", rst2 - rst1);
  fprintf(stderr, "Elapsed time(s) : %ld\n", t2 - t1);

  return(0);
}

///////////////////////////////////////
// Function: trim - Remove "\n"      //
///////////////////////////////////////

int trim(char *line) {
  int end_pos = strlen(line) - 1;

  if (line[end_pos] == '\n') {
    line[end_pos] = '\0';
    return 1;
  }
  return 0;
}

///////////////////////////////////////////////////////
// Function: get_genome_inf - Get genome information //
///////////////////////////////////////////////////////

int get_genome_inf() {
  FILE *fp, *fp_out;
  char line[BUF_SIZE];
  int ret;
  long max_len = 0;

  fprintf(stderr, ":::: Reference stats ::::\n\n");
  fprintf(stderr, "file name : %s\n", genome.file);
  fprintf(stderr, "\n");

  if ((fp = fopen(genome.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", genome.file);
    return FAILED;
  }

  genome.num_seq = 0;
  genome.len = 0;

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      if (genome.num_seq != 0) {
        if (genome.len < REF_SEQ_LEN_MIN) {
          fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %d.\n", REF_SEQ_LEN_MIN);
          return FAILED;
        }
        fprintf(stderr, "ref.%ld (len:%ld) : %s\n", genome.num_seq, genome.len, genome.id);
        fclose(fp_out);
        if (genome.len > max_len) {
          max_len = genome.len;
        }
      }

      genome.num_seq ++;
      if (genome.num_seq > REF_SEQ_NUM_MAX) {
        fprintf(stderr, "ERROR: References are too many. Max number of reference is %d.\n", REF_SEQ_NUM_MAX);
        return FAILED;
      }

      strncpy(genome.id, line + 1, REF_ID_LEN_MAX);
      genome.id[REF_ID_LEN_MAX] = '\0';

      sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, genome.num_seq);
      if ((fp_out = fopen(sim.outfile_ref, "w")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open output file: %s\n", sim.outfile_ref);
        return FAILED;
      }

      genome.len = 0;

      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }

      fprintf(fp_out, ">%s\n", genome.id);
    } else {
      genome.len += strlen(line);

      if (genome.len > REF_SEQ_LEN_MAX) {
        fprintf(stderr, "ERROR: Reference is too long. Acceptable length <= %d.\n", REF_SEQ_LEN_MAX);
        return FAILED;
      }

      fprintf(fp_out, "%s\n", line);
    }
  }
  fclose(fp);

  if (genome.len < REF_SEQ_LEN_MIN) {
    fprintf(stderr, "ERROR: Reference is too short. Acceptable length >= %d.\n", REF_SEQ_LEN_MIN);
    return FAILED;
  }
  fprintf(stderr, "ref.%ld (len:%ld) : %s\n", genome.num_seq, genome.len, genome.id);
  fclose(fp_out);
  if (genome.len > max_len) {
    max_len = genome.len;
  }

  fprintf(stderr, "\n");

  if ((genome.seq = (char *)malloc(max_len + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((genome.hp = (short *)malloc(max_len * sizeof(short) + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////
// Function: get_genome_seq - Get genome sequence //
////////////////////////////////////////////////////

int get_genome_seq() {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;
  long i, j;
  long nstart, nend;
  short nnum;

  sprintf(sim.outfile_ref, "%s_%04ld.ref", sim.prefix, genome.num);

  if ((fp = fopen(sim.outfile_ref, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", sim.outfile_ref);
    return FAILED;
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(genome.seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  genome.seq[offset] = '\0';
  genome.len = strlen(genome.seq);

  for(i=0; i<genome.len; i++) {
    genome.seq[i] = toupper(genome.seq[i]);
  }

  nstart = 0;
  nend = 0;
  nnum = 1; 
  for(i=1; i<=genome.len; i++) {
    if ((i < genome.len ) && (genome.seq[i-1] == genome.seq[i])) {
      nend = i;
      nnum ++;
      if (nnum > 11) {
        nnum = 10;
      }
    } else {
      if (genome.seq[i-1] == 'N') {
        for(j=nstart; j<=nend; j++) {
          genome.hp[j] = 1;
          genome.hpfreq[1] ++;
        }
      } else {
        for(j=nstart; j<=nend; j++) {
          genome.hp[j] = nnum;
          genome.hpfreq[nnum] ++;
        }
      }
      nstart = i;
      nend = nstart;
      nnum = 1;
    }
  }

  return SUCCEEDED;
}


///////////////////////////////////////////////////////////////
// Function: get_transcript_inf - Get transcript information //
///////////////////////////////////////////////////////////////

int get_transcript_inf() {
  FILE *fp;
  char *tp;
  char line[BUF_SIZE];
  int flg1, flg2;
  int ret;
  long max_len = 0;

  transcript.num_seq = 0;
  transcript.total_exp = 0.0;
  transcript.rank_max = 0;

  if ((fp = fopen(transcript.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", transcript.file);
    return FAILED;
  }

  transcript.len = 0;
  flg1 = 1;

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    if (trim(line) == EXISTS_LINE_FEED) {
      flg2 = 1;
    } else {
      flg2 = 0;
    }
    if (flg1 == 1) {
      tp = strtok(line, "\t");
      transcript.num_seq ++;
      tp = strtok(NULL, "\t");
      transcript.total_exp += atoi(tp);
      tp = strtok(NULL, "\t");
      transcript.total_exp += atoi(tp);
      tp = strtok(NULL, "\t");
      transcript.len = strlen(tp);

    } else {
      transcript.len += strlen(line);
    }
    if (flg2 == 1) {
      if (transcript.len > max_len) {
        max_len = transcript.len;
      }
    }
    flg1 = flg2;
  }
  fclose(fp);

  if ((transcript.seq = (char *)malloc(max_len + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((transcript.hp = (short *)malloc(max_len * sizeof(short) + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  transcript.rank_max = ceil((float)max_len / 1000);

  return SUCCEEDED;
}

///////////////////////////////////////////////////////////////
// Function: print_transcript_stats - Print transcript stats //
///////////////////////////////////////////////////////////////

void print_transcript_stats() {
  fprintf(stderr, ":::: transcript stats ::::\n\n");

  fprintf(stderr, "file name : %s\n", transcript.file);
  fprintf(stderr, "transcript num : %ld\n", transcript.num_seq);
  fprintf(stderr, "total expression value : %ld\n", transcript.total_exp);
  fprintf(stderr, "\n");
}

//////////////////////////////////////////////////////
// Function: get_sample_inf - Get FASTQ information //
//////////////////////////////////////////////////////

int get_sample_inf() {
  FILE *fp;
  char *tp, *item;
  char line[BUF_SIZE];
  char qc_tmp[FASTQ_LEN_MAX];
  long len;
  double prob;
  double accuracy;
  double accuracy_total = 0;
  long value;
  double variance;
  long i;
  int line_num;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  sample.num = 0;
  sample.len_min = LONG_MAX;
  sample.len_max = 0;
  sample.len_total = 0;
  sample.num_filtered = 0;
  sample.len_min_filtered = LONG_MAX;
  sample.len_max_filtered = 0;
  sample.len_total_filtered = 0;

  if (sim.method == METHOD_SAM_REUSE) {
    while (fgets(line, BUF_SIZE, fp_stats) != NULL) {
      trim(line);
      tp = strtok(line, "\t");
      item = tp;
      tp = strtok(NULL, "\t");

      if (strcmp(item, "num") == 0) {
        sample.num_filtered = atol(tp);
      } else if (strcmp(item, "len_total") == 0) {
        sample.len_total_filtered = atol(tp);
      } else if (strcmp(item, "len_min") == 0) {
        sample.len_min_filtered = atol(tp);
      } else if (strcmp(item, "len_max") == 0) {
        sample.len_max_filtered = atol(tp);
      } else if (strcmp(item, "len_mean") == 0) {
        sample.len_mean_filtered = atof(tp);
      } else if (strcmp(item, "len_sd") == 0) {
        sample.len_sd_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_mean") == 0) {
        sample.accuracy_mean_filtered = atof(tp);
      } else if (strcmp(item, "accuracy_sd") == 0) {
        sample.accuracy_sd_filtered = atof(tp);
      }
    }
  } else {
    if ((fp = fopen(sample.file, "r")) == NULL) {
      fprintf(stderr, "ERROR: Cannot open file: %s\n", sample.file);
      return FAILED;
    }

    qc_tmp[0] = '\0';
    len = 0;
    line_num = 0;

    while (fgets(line, BUF_SIZE, fp) != NULL) {
      if (trim(line) == EXISTS_LINE_FEED) {
        line_num ++;

        if (line_num == 4) {
          len += strlen(line);

          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %d.\n", FASTQ_LEN_MAX);
            return FAILED;
          }

          sample.num ++;
          sample.len_total += len;

          if (sample.num > FASTQ_NUM_MAX) {
            fprintf(stderr, "ERROR: fastq is too many. Max acceptable number is %d.\n", FASTQ_NUM_MAX);
            return FAILED;
          }

          if (len > sample.len_max) {
            sample.len_max = len;
          }
          if (len < sample.len_min) {
            sample.len_min = len;
          }

          if ((len >= sim.len_min) && (len <= sim.len_max)) {
            strcat(qc_tmp, line);
            prob = 0.0;
            for (i=0; i<len; i++) {
              prob += qc[(int)qc_tmp[i] - 33].prob;
            }
            accuracy = 1.0 - (prob / len);

            if ((accuracy >= sim.accuracy_min) && (accuracy <= sim.accuracy_max)) {
              accuracy_total += accuracy;
              sample.num_filtered ++;
              sample.len_total_filtered += len;

              freq_len[len] ++;
              value = (int)(accuracy * 100000 + 0.5); 
              freq_accuracy[value] ++;

              fprintf(fp_filtered, "%s\n", qc_tmp);

              if (len > sample.len_max_filtered) {
                sample.len_max_filtered = len;
              }
              if (len < sample.len_min_filtered) {
                sample.len_min_filtered = len;
              }
            }
          }

          line_num = 0;
          qc_tmp[0] = '\0';
          len = 0;
        }
      } else {
        if (line_num == 3) {
          len += strlen(line);
          if (len > FASTQ_LEN_MAX) {
            fprintf(stderr, "ERROR: fastq is too long. Max acceptable length is %d.\n", FASTQ_LEN_MAX);
            return FAILED;
          }
          strcat(qc_tmp, line);
        }
      }
    }

    fclose(fp);

    if (sample.num_filtered < 1) {
      fprintf(stderr, "ERROR: there is no sample in the valid range of length and accuracy.\n");
      return FAILED;
    }

    sample.len_mean_filtered = (double)sample.len_total_filtered / sample.num_filtered;
    sample.accuracy_mean_filtered = accuracy_total / sample.num_filtered;

    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) { 
        variance += pow((sample.len_mean_filtered - i), 2) * freq_len[i];
      }
    }
    sample.len_sd_filtered = sqrt(variance / sample.num_filtered);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) { 
        variance += pow((sample.accuracy_mean_filtered - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sample.accuracy_sd_filtered = sqrt(variance / sample.num_filtered);

    if (sim.method == METHOD_SAM_STORE) {
      fprintf(fp_stats, "num\t%ld\n", sample.num_filtered);
      fprintf(fp_stats, "len_total\t%lld\n", sample.len_total_filtered);
      fprintf(fp_stats, "len_min\t%ld\n", sample.len_min_filtered);
      fprintf(fp_stats, "len_max\t%ld\n", sample.len_max_filtered);
      fprintf(fp_stats, "len_mean\t%f\n", sample.len_mean_filtered);
      fprintf(fp_stats, "len_sd\t%f\n", sample.len_sd_filtered);
      fprintf(fp_stats, "accuracy_mean\t%f\n", sample.accuracy_mean_filtered);
      fprintf(fp_stats, "accuracy_sd\t%f\n", sample.accuracy_sd_filtered);
    }
  }

  return SUCCEEDED;
}

//////////////////////////////////////////////////////
// Function: print_sample_stats - Print FASTQ stats //
//////////////////////////////////////////////////////

void print_sample_stats() {
  fprintf(stderr, ":::: sample reads stats ::::\n\n");

  if (sim.method == METHOD_SAM_REUSE) {
    fprintf(stderr, "file name : %s\n", sim.profile_fq);
  } else {
    fprintf(stderr, "file name : %s\n", sample.file);
    fprintf(stderr, "\n:: all reads ::\n");
    fprintf(stderr, "read num. : %ld\n", sample.num);
    fprintf(stderr, "read total length : %lld\n", sample.len_total);
    fprintf(stderr, "read min length : %ld\n", sample.len_min);
    fprintf(stderr, "read max length : %ld\n", sample.len_max);
  }

  fprintf(stderr, "\n:: filtered reads ::\n");
  fprintf(stderr, "read num. : %ld\n", sample.num_filtered);
  fprintf(stderr, "read total length : %lld\n", sample.len_total_filtered);
  fprintf(stderr, "read min length : %ld\n", sample.len_min_filtered);
  fprintf(stderr, "read max length : %ld\n", sample.len_max_filtered);
  fprintf(stderr, "read length mean (SD) : %f (%f)\n",
    sample.len_mean_filtered, sample.len_sd_filtered);
  fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
    sample.accuracy_mean_filtered, sample.accuracy_sd_filtered);
  fprintf(stderr, "\n");
}

////////////////////////////////////////////////////////
// Function: get_templ_inf - Get Template information //
////////////////////////////////////////////////////////

int get_templ_inf() {
  FILE *fp;
  char line[BUF_SIZE];
  long seqlen, len;
  int ret;

  templ.num = 0;
  templ.len_total = 0;

  if ((fp = fopen(templ.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", templ.file);
    return FAILED;
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    ret = trim(line);

    if (line[0] == '>') {
      templ.num ++;
      if (templ.num > TEMPLATE_NUM_MAX) {
        fprintf(stderr, "ERROR: template is too many. Max acceptable number is %d.\n", TEMPLATE_NUM_MAX);
        return FAILED;
      }
      seqlen = 0;
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      len = strlen(line);
      seqlen += len;
      templ.len_total += len;
      if (seqlen > TEMPLATE_LEN_MAX) {
        fprintf(stderr, "ERROR: template is too long. Max acceptable length is %d.\n", TEMPLATE_LEN_MAX);
        return FAILED;
      }
    }
  }

  if ((templ.seq = (char *)malloc(TEMPLATE_LEN_MAX + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((templ.hp = (short *)malloc(TEMPLATE_LEN_MAX * sizeof(short) + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////
// Function: print_templ_stats - Print Template stats //
////////////////////////////////////////////////////////

void print_templ_stats() {
  fprintf(stderr, ":::: Template stats ::::\n\n");

  fprintf(stderr, "file name : %s\n", templ.file);
  fprintf(stderr, "template num. : %ld\n", templ.num);
  fprintf(stderr, "template total length : %lld\n", templ.len_total);
  fprintf(stderr, "\n");
}

//////////////////////////////////////////////////////////
// Function: init_sim_res - Initiate simulation results //
//////////////////////////////////////////////////////////

void init_sim_res() {
  sim.res_num = 0;
  sim.res_len_total = 0;
  sim.res_sub_num = 0;
  sim.res_ins_num = 0;
  sim.res_del_num = 0;
  sim.res_len_min = LONG_MAX;
  sim.res_len_max = 0;
}

/////////////////////////////////////////////////////////
// Function: set_sim_param - Set simulation parameters //
/////////////////////////////////////////////////////////

int set_sim_param() {
  FILE *fp;
  long sum;

  // strategy and method
  if ((!(sim.set_flg[0])) || (!(sim.set_flg[1]))) {
    fprintf(stderr, "ERROR: --strategy and --method must be set.\n");
    return FAILED;
  }

  if ((sim.strategy != STRATEGY_WGS) && (sim.method == METHOD_SAM)) {
    fprintf(stderr, "ERROR: sampling-based simulation is possible only for wgs strategy.\n");
    return FAILED;
  }

  // WGS
  if (sim.strategy == STRATEGY_WGS) {
    if (!(sim.set_flg[2])) {
      fprintf(stderr, "ERROR: for --strategy wgs, --genome must be set.\n");
      return FAILED;
    }
  // transcript
  } else if (sim.strategy == STRATEGY_TRANS) {
    if (!(sim.set_flg[3])) {
      fprintf(stderr, "ERROR: for --strategy trans, --transcript must be set.\n");
      return FAILED;
    }
  // template
  } else {
    if (!(sim.set_flg[21])) {
      fprintf(stderr, "ERROR: for --strategy templ, --template must be set.\n");
      return FAILED;
    }
  }

  // prefix and outfile
  if (!(sim.set_flg[4])) {
    if ((sim.prefix = (char *)malloc(3)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(sim.prefix, "sd");
  }

  if ((sim.outfile_ref = (char *)malloc(strlen(sim.prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim.outfile_fq = (char *)malloc(strlen(sim.prefix) + 12)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim.outfile_maf = (char *)malloc(strlen(sim.prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((sim.outfile_sam = (char *)malloc(strlen(sim.prefix) + 10)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  // prefix of read ID
  if (!(sim.set_flg[5])) {
    if ((sim.id_prefix = (char *)malloc(2)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      exit(-1);
    }
    strcpy(sim.id_prefix, "S");
  }

  // depth
  if (!(sim.set_flg[6])) {
    sim.depth = 20.0;
  }

  // length-min
  if (!(sim.set_flg[7])) {
    sim.len_min = 100;
  }

  // length-max
  if (!(sim.set_flg[8])) {
    sim.len_max = 1000000;
  }

  // difference-ratio
  if (!(sim.set_flg[9])) {
    sim.sub_ratio = 6;
    sim.ins_ratio = 55;
    sim.del_ratio = 39;
  }

  sum = sim.sub_ratio + sim.ins_ratio + sim.del_ratio;
  sim.sub_rate = (double)sim.sub_ratio / sum;
  sim.ins_rate = (double)sim.ins_ratio / sum;
  sim.del_rate = (double)sim.del_ratio / sum;

  //sample, sample-profile-id
  if (sim.method == METHOD_SAM) {
    if (sim.set_flg[11]) {
      if (sim.set_flg[12]) {
        sim.method = METHOD_SAM_STORE;
      }
    } else {
      if (sim.set_flg[12]) {
        sim.method = METHOD_SAM_REUSE;
      } else {
        fprintf(stderr, "ERROR: for --method sample, --sample (and/or --sample-profile-id) must be set.\n");
        return FAILED;
      }
    }
  }
  if (sim.set_flg[12]) {
    if ((sim.profile_fq = (char *)malloc(strlen(sim.profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }
    if ((sim.profile_stats = (char *)malloc(strlen(sim.profile_id) + 22)) == 0) {
      fprintf(stderr, "ERROR: Cannot allocate memory.\n");
      return FAILED;
    }
    sprintf(sim.profile_fq, "sample_profile_%s.fastq", sim.profile_id);
    sprintf(sim.profile_stats, "sample_profile_%s.stats", sim.profile_id);
  }
  if (sim.method == METHOD_SAM_STORE) {
    if ((fp = fopen(sim.profile_fq, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim.profile_fq);
      fclose(fp);
      return FAILED;
    }
    if ((fp = fopen(sim.profile_stats, "r")) != NULL) {
      fprintf(stderr, "ERROR: %s exists.\n", sim.profile_stats);
      fclose(fp);
      return FAILED;
    }
  }
  if (sim.method == METHOD_SAM_REUSE) {
    if ((fp = fopen(sim.profile_fq, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim.profile_fq);
      return FAILED;
    }
    fclose(fp);
    if ((fp = fopen(sim.profile_stats, "r")) == NULL) {
      fprintf(stderr, "ERROR: %s does not exist.\n", sim.profile_stats);
      return FAILED;
    }
    fclose(fp);
  }

  // accuracy-min
  if (sim.set_flg[13]) {
    sim.accuracy_min = int(sim.accuracy_min * 100) * 0.01;
  } else {
    sim.accuracy_min = 0.75;
  }

  // accuracy-max
  if (sim.set_flg[14]) {
    sim.accuracy_max = int(sim.accuracy_max * 100) * 0.01;
  } else {
    sim.accuracy_max = 1.0;
  }

  // qshmm
  if (sim.method == METHOD_QS) {
    if ((!sim.set_flg[15])) {
      fprintf(stderr, "ERROR: for --method qshmm, --qshmm must be set.\n");
      return FAILED;
    }
  }

  // errhmm
  if (sim.method == METHOD_ERR) {
    if ((!sim.set_flg[16])) {
      fprintf(stderr, "ERROR: for --method errhmm, --errhmm must be set.\n");
      return FAILED;
    }
  }

  // length-mean
  if (!(sim.set_flg[17])) {
    sim.len_mean = 9000;
  }

  // length-sd
  if (!(sim.set_flg[18])) {
    sim.len_sd = 7000;
  }

  // accuracy-mean
  if (sim.set_flg[19]) {
    sim.accuracy_mean = int(sim.accuracy_mean * 100) * 0.01;
  } else {
    sim.accuracy_mean = 0.85;
  }

  // length and accuracy
  if (sim.len_min > sim.len_max) {
    fprintf(stderr, "ERROR: length min(%ld) is greater than max(%ld).\n", sim.len_min, sim.len_max);
    return FAILED;
  }

  // pass-num
  if (!(sim.set_flg[20])) {
    sim.pass_num = 1;
  }
  if (sim.pass_num > 1) {
    if ((sim.method == METHOD_SAM) || (sim.method == METHOD_SAM_STORE) || (sim.method == METHOD_SAM_REUSE)) {
      fprintf(stderr, "ERROR: sampling-based simulation supports only single-pass.\n");
      return FAILED;
    }
  }

  // hp-del-bias
  if (!(sim.set_flg[22])) {
    sim.hp_del_bias = 1;
  }

  return SUCCEEDED;
}

///////////////////////////////////////////////////////
// Function: simulate_by_sample - Simulate by sample //
///////////////////////////////////////////////////////

int simulate_by_sample() {
  long len;
  long long len_total = 0;
  long sample_num, sample_interval, sample_value, sample_residue;
  long num;
  long i, j;
  long index;
  long value;
  double accuracy, accuracy_total = 0.0;
  double prob, variance;
  char id[128], nt;
  int digit_num1[4], digit_num2[4], digit_num[4];
  long read_offset, ref_offset, maf_offset;
  long qc_value, rand_value;
  int hp;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  sample_num = (long)(sim.len_quota / sample.len_total_filtered);
  sample_residue = sim.len_quota % sample.len_total_filtered;
  if (sample_residue == 0) {
    sample_interval = 1;
  } else {
    sample_interval = (long)((double)(sample.len_total_filtered / sample_residue) * 2 + 0.5);
    if (sample_interval > (long)(sample.num_filtered * 0.5)) {
      sample_interval = (long)(sample.num_filtered * 0.5);
    }
  }

  // Make simulation data
  while (len_total < sim.len_quota) {
    rewind(fp_filtered);

    sample_value = rand() % sample.num_filtered;
    while (fgets(mut.qc, sample.len_max_filtered + 2, fp_filtered) != NULL) {
      if (len_total >= sim.len_quota) {
        break;
      }

      trim(mut.qc);

      if (sample_value % sample_interval == 0) {
        num = sample_num + 1;
      } else {
        num = sample_num;
      }
      sample_value ++;

      for (i=0; i<num; i++) {
        if (len_total >= sim.len_quota) {
          break;
        }

        mut.len = strlen(mut.qc);

        if (mut.len >= genome.len) {
          mut.offset = 0;
          mut.len = genome.len;
        } else {
          mut.offset = rand() % (genome.len - mut.len + 1);
        }

        sim.res_num ++;

        for (j=0; j<mut.len; j++) {
          mut.seq[j] = genome.seq[mut.offset + j];
          mut.hp[j] = genome.hp[mut.offset + j];
        }
        mut.seq[mut.len] = '\0';
        if (sim.res_num % 2 == 1) {
          mut.seq_strand = '+';
        } else {
          mut.seq_strand = '-';
          revcomp(mut.seq);
          revshort(mut.hp, mut.len);
        }

        ref_offset = 0;
        read_offset = 0;
        maf_offset = 0;
        while ((ref_offset < mut.len) && (read_offset < mut.len)) {
          nt = mut.seq[ref_offset];
          qc_value = (int)mut.qc[read_offset] - 33;
          rand_value = rand() % 1000000;
          if (rand_value < mut.sub_thre[qc_value]) {
            sim.res_sub_num ++;
            index = rand() % 3;
            if (nt == 'A') {
              mut.read_seq[read_offset] = mut.sub_nt_a[index];
            } else if (nt == 'T') {
              mut.read_seq[read_offset] = mut.sub_nt_t[index];
            } else if (nt == 'G') {
              mut.read_seq[read_offset] = mut.sub_nt_g[index];
            } else if (nt == 'C') {
              mut.read_seq[read_offset] = mut.sub_nt_c[index];
            } else {
              index = rand() % 4;
              mut.read_seq[read_offset] = mut.sub_nt_n[index];
            }
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
          } else if (rand_value < mut.ins_thre[qc_value]) {
            sim.res_ins_num ++;
            index = rand() % 8;
            if (index >= 4) {
              mut.read_seq[read_offset] = nt;
            } else {
              mut.read_seq[read_offset] = mut.ins_nt[index];
            }
            mut.maf_ref_seq[maf_offset] = '-';
          } else {
            mut.read_seq[read_offset] = nt;
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
          }
          mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
          maf_offset ++;
          read_offset ++;

          while ((ref_offset < mut.len) && (read_offset < mut.len)) {
            hp = mut.hp[ref_offset-1];
            rand_value = rand() % 1000000;
            qc_value = (int)mut.qc[read_offset-1] - 33;
            if (rand_value < mut.del_thre[qc_value] * genome.hp_del_bias[hp]) {
              sim.res_del_num ++;
              mut.maf_seq[maf_offset] = '-';
              mut.maf_ref_seq[maf_offset] = mut.seq[ref_offset];
              maf_offset ++;
              ref_offset ++;
            } else {
              break;
            }
          }
        }
        mut.qc[read_offset] = '\0';
        mut.read_seq[read_offset] = '\0';
        mut.maf_seq[maf_offset] = '\0';
        mut.maf_ref_seq[maf_offset] = '\0';

        if (mut.seq_strand == '-') {
          revcomp(mut.maf_seq);
          revcomp(mut.maf_ref_seq);
        }

        mut.seq_left = mut.offset + 1;
        mut.seq_right = mut.offset + ref_offset;

        len = strlen(mut.read_seq);
        sim.res_len_total += len;
        len_total += len;
        freq_len[len] ++;

        if (len > sim.res_len_max) {
          sim.res_len_max = len;
        }
        if (len < sim.res_len_min) {
          sim.res_len_min = len;
        }

        prob = 0.0;
        for (j=0; j<len; j++) {
          prob += qc[(int)mut.qc[j] - 33].prob;
        }
        accuracy = 1.0 - (prob / len);
        accuracy_total += accuracy;
        value = (int)(accuracy * 100000 + 0.5);
        freq_accuracy[value] ++;

        sprintf(id, "%s%ld_%ld", sim.id_prefix, genome.num, sim.res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.qc);

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim.res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut.seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(genome.len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns ref"); 
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut.seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", genome.len, mut.maf_ref_seq);
        fprintf(fp_maf, "s %s", id); 
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
      }
    }

    sample_num = 0;
  }

  sim.res_len_mean = (double)sim.res_len_total / sim.res_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_num;

  if (sim.res_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////
// Function: simulate_by_qshmm - Simulate by Model //
/////////////////////////////////////////////////////

int simulate_by_qshmm() {
  long len;
  long long len_total = 0;
  long h, i, j, k, l;
  long state;
  double prob, mean, variance, sd;
  double kappa, theta, gamma;
  double len_prob_total, freq_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  static long prob2len[100001], prob2accuracy[100001];
  static long freq2qc[ACCURACY_MAX+1][1001];
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2qc[ACCURACY_MAX+1][STATE_MAX+1][101];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][101];
  long len_rand_value, accuracy_rand_value;
  long qc_rand_value_freq[ACCURACY_MAX+1];
  long qc_rand_value_init[ACCURACY_MAX+1];
  long qc_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long qc_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, pre_index;
  long acc_wk, accuracy_min, accuracy_max;
  char id[128], nt;
  int digit_num1[4], digit_num2[4], digit_num[4];
  int qeval;
  long read_offset, ref_offset, maf_offset;
  long qc_value, rand_value;
  int hp;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // length distribution
  variance = pow(sim.len_sd, 2);
  kappa = pow(sim.len_mean, 2) / variance; 
  theta = variance / sim.len_mean;
  gamma = tgamma(kappa);

  if (sim.len_sd == 0.0) {
    prob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1; 
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += pow(i, kappa-1) * exp(-1 * i / theta) / pow(theta, kappa) / gamma;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (sim.pass_num == 1) {
    if (len_rand_value < 1) {
      fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
      return FAILED;
    }
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (qshmm.exist_hmm[i] == 1) {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=1; j<=STATE_MAX; j++) {
        if (qshmm.ip[i][j] == 0) {
          continue;
        }
        qc_prob_total += qshmm.ip[i][j];
        end_wk = int(qc_prob_total * 100 + 0.5);
        if (end_wk > 100) {
          end_wk = 100;
        }

        for (k=start_wk; k<=end_wk; k++) {
          init2state[i][k] = j;
        }

        if (end_wk >= 100) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_init[i] = end_wk;

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=0; k<=93; k++) {
          if (qshmm.ep[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.ep[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            emis2qc[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_emis[i][j] = end_wk;
      }

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=1; k<=STATE_MAX; k++) {
          if (qshmm.tp[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.tp[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            tran2state[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_tran[i][j] = end_wk;
      }
    } else {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=0; j<=93; j++) {
        if (uni_ep[i][j] == 0) {
          continue;
        }
        qc_prob_total += uni_ep[i][j];
        end_wk = int(qc_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (k=start_wk; k<=end_wk; k++) {
          freq2qc[i][k] = j;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_freq[i] = end_wk;
    }
  }

  // simulation
  while (len_total < sim.len_quota) {
    index = rand() % len_rand_value + 1;
    mut.len = prob2len[index];
    if (len_total + mut.len > sim.len_quota) {
      mut.len = sim.len_quota - len_total;
      if (mut.len < sim.len_min) {
        mut.len = sim.len_min;
      }
    }

    index = rand() % accuracy_rand_value + 1;
    mut.acc = prob2accuracy[index];
    if (mut.len >= genome.len) {
      mut.offset = 0;
      mut.len = genome.len;
    } else {
      mut.offset = rand() % (genome.len - mut.len + 1);
    }

    mut.seq_left = mut.offset + 1;
    mut.seq_right = mut.offset + mut.len;
    sim.res_num ++;

    for (i=0; i<mut.len; i++) {
      mut.seq[i] = genome.seq[mut.offset + i];
      mut.hp[i] = genome.hp[mut.offset + i];
    }
    mut.seq[mut.len] = '\0';
    if (sim.res_num % 2 == 1) {
      mut.seq_strand = '+';
    } else {
      mut.seq_strand = '-';
      revcomp(mut.seq);
      revshort(mut.hp, mut.len);
    }

    for (h=0; h<sim.pass_num; h++) {
      ref_offset = 0;
      read_offset = 0;
      maf_offset = 0;
      while (ref_offset < mut.len) {
        if (qshmm.exist_hmm[mut.acc] == 1) {
          if (read_offset == 0) {
            index = rand() % qc_rand_value_init[mut.acc] + 1;
            state = init2state[mut.acc][index];
          } else {
            index = rand() % qc_rand_value_tran[mut.acc][state] + 1;
            state = tran2state[mut.acc][state][index];
          }
          index = rand() % qc_rand_value_emis[mut.acc][state] + 1;
          index = emis2qc[mut.acc][state][index];
        } else {
          index = rand() % qc_rand_value_freq[mut.acc] + 1;
          index = freq2qc[mut.acc][index];
        }
        mut.qc[read_offset] = qc[index].character;

        nt = mut.seq[ref_offset];
        qc_value = (int)mut.qc[read_offset] - 33;
        rand_value = rand() % 1000000;
        if (rand_value < mut.sub_thre[qc_value]) {
          sim.res_sub_num ++;
          index = rand() % 3;
          if (nt == 'A') {
            mut.read_seq[read_offset] = mut.sub_nt_a[index];
          } else if (nt == 'T') {
            mut.read_seq[read_offset] = mut.sub_nt_t[index];
          } else if (nt == 'G') {
            mut.read_seq[read_offset] = mut.sub_nt_g[index];
          } else if (nt == 'C') {
            mut.read_seq[read_offset] = mut.sub_nt_c[index];
          } else {
            index = rand() % 4;
            mut.read_seq[read_offset] = mut.sub_nt_n[index];
          }
          mut.maf_ref_seq[maf_offset] = nt;
          ref_offset ++;
        } else if (rand_value < mut.ins_thre[qc_value]) {
          sim.res_ins_num ++;
          index = rand() % 8;
          if (index >= 4) {
            mut.read_seq[read_offset] = nt;
          } else {
            mut.read_seq[read_offset] = mut.ins_nt[index];
          }
          mut.maf_ref_seq[maf_offset] = '-';
        } else {
          mut.read_seq[read_offset] = nt;
          mut.maf_ref_seq[maf_offset] = nt;
          ref_offset ++;
        }
        mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
        maf_offset ++;
        read_offset ++;

        while (ref_offset < mut.len) {
          hp = mut.hp[ref_offset-1];
          rand_value = rand() % 1000000;
          qc_value = (int)mut.qc[read_offset-1] - 33;
          if (rand_value < mut.del_thre[qc_value] * genome.hp_del_bias[hp]) {
            sim.res_del_num ++;
            mut.maf_seq[maf_offset] = '-';
            mut.maf_ref_seq[maf_offset] = mut.seq[ref_offset];
            maf_offset ++;
            ref_offset ++;
          } else {
            break;
          }
        }
      }
      mut.qc[read_offset] = '\0';
      mut.read_seq[read_offset] = '\0';
      mut.maf_seq[maf_offset] = '\0';
      mut.maf_ref_seq[maf_offset] = '\0';

      if (mut.seq_strand == '-') {
        revcomp(mut.maf_seq);
        revcomp(mut.maf_ref_seq);
      }

      len = strlen(mut.read_seq);
      sim.res_len_total += len;

      if (h == 0) {
        len_total += len;
      }

      freq_len[len] ++;

      if (len > sim.res_len_max) {
        sim.res_len_max = len;
      }
      if (len < sim.res_len_min) {
        sim.res_len_min = len;
      }

      prob = 0.0;
      for (i=0; i<len; i++) {
        prob += qc[(int)mut.qc[i] - 33].prob;
      }
      value = 1.0 - (prob / len);
      accuracy_total += value;
      acc_wk = (int)(value * 100000 + 0.5);
      freq_accuracy[acc_wk] ++;

      if (sim.pass_num == 1) {
        sprintf(id, "%s%ld_%ld", sim.id_prefix, genome.num, sim.res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.qc);
      } else {
        sprintf(id, "%s%ld/%ld/%ld", sim.id_prefix, genome.num, sim.res_num, h);
        fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.qc);
        fprintf(fp_sam, "\tcx:i:3\tip:B:C");
        for (i=0; i<len; i++) {
          fprintf(fp_sam, ",9");
        }
        fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
        for (i=0; i<len; i++) {
          fprintf(fp_sam, ",9");
        }
        qeval = len - 1;
        fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
      }

      digit_num1[0] = 3;
      digit_num2[0] = 1 + count_digit(sim.res_num);
      digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

      digit_num1[1] = count_digit((mut.seq_left - 1));
      digit_num2[1] = 1;
      digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

      digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
      digit_num2[2] = count_digit(len);
      digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

      digit_num1[3] = count_digit(genome.len);
      digit_num2[3] = count_digit(len);
      digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

      fprintf(fp_maf, "a\ns ref");
      while (digit_num1[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num1[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld", mut.seq_left - 1);
      while (digit_num1[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
      while (digit_num1[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n", genome.len, mut.maf_ref_seq);
      fprintf(fp_maf, "s %s", id);
      while (digit_num2[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num2[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %d", 0);
      while (digit_num2[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
      while (digit_num2[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
    }
  }

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}

//////////////////////////////////////////////////////////////////////////
// Function: simulate_by_qshmm_trans - Simulate by qshmm for transcript //
//////////////////////////////////////////////////////////////////////////

int simulate_by_qshmm_trans() {
  FILE *fp;
  char *tp;
  char line[BUF_SIZE];
  int flg1, flg2;
  long offset, copy_size;
  long len;
  long h, i, j, k, l;
  double prob, mean, variance, sd;
  double kappa, theta, gamma;
  double len_prob_total, freq_total, accuracy_prob_total, qc_prob_total;
  double ssp_prob_total;
  double value, sum;
  double accuracy_total = 0.0;
  static long prob2len[100001], prob2accuracy[100001];
  static long prob2ssp[TR_RANK_MAX][1001];
  static long freq2qc[ACCURACY_MAX+1][1001];
  long state;
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2qc[ACCURACY_MAX+1][STATE_MAX+1][101];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][101];
  long len_rand_value, accuracy_rand_value;
  long ssp_rand_value[TR_RANK_MAX+1];
  long qc_rand_value_freq[ACCURACY_MAX+1];
  long qc_rand_value_init[ACCURACY_MAX+1];
  long qc_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long qc_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, pre_index;
  long acc_wk, accuracy_min, accuracy_max;
  char id[128], nt;
  int digit_num1[4], digit_num2[4], digit_num[4];
  int qeval;
  int read_num, rank;
  long read_offset, ref_offset, maf_offset;
  long qc_value, rand_value;
  int hp;
  long nstart, nend;
  short nnum;
  long sum1, sum2;
  double rate;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // length distribution
  variance = pow(sim.len_sd, 2);
  kappa = pow(sim.len_mean, 2) / variance;
  theta = variance / sim.len_mean;
  gamma = tgamma(kappa);

  if (sim.len_sd == 0.0) {
    prob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1;
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += pow(i, kappa-1) * exp(-1 * i / theta) / pow(theta, kappa) / gamma;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (len_rand_value < 1) {
    fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
    return FAILED;
  }

  // sequencing start pos distribution
  for (i=1; i<=transcript.rank_max; i++) {
    sum = 0;
    value = (double)1 / i;
    for (j=1; j<=21; j++) {
      sum += value / pow(j,(1+value));
    }
    start_wk = 1;
    ssp_prob_total = 0.0;
    for (j=1; j<=21; j++) {
      ssp_prob_total += (value / pow(j,(1+value))) / sum;
      end_wk = int(ssp_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }
      for (k=start_wk; k<=end_wk; k++) {
        prob2ssp[i][k] = (j-1) * 5;
      }
      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    ssp_rand_value[i] = end_wk;
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1;
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }
    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }
    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (qshmm.exist_hmm[i] == 1) {
      start_wk = 1;
      qc_prob_total = 0.0;

      for (j=1; j<=STATE_MAX; j++) {
        if (qshmm.ip[i][j] == 0) {
          continue;
        }
        qc_prob_total += qshmm.ip[i][j];
        end_wk = int(qc_prob_total * 100 + 0.5);
        if (end_wk > 100) {
          end_wk = 100;
        }
        for (k=start_wk; k<=end_wk; k++) {
          init2state[i][k] = j;
        }
        if (end_wk >= 100) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_init[i] = end_wk;

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1;
        qc_prob_total = 0.0;

        for (k=0; k<=93; k++) {
          if (qshmm.ep[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.ep[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }
          for (l=start_wk; l<=end_wk; l++) {
            emis2qc[i][j][l] = k;
          }
          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_emis[i][j] = end_wk;
      }

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1;
        qc_prob_total = 0.0;

        for (k=1; k<=STATE_MAX; k++) {
          if (qshmm.tp[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.tp[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }
          for (l=start_wk; l<=end_wk; l++) {
            tran2state[i][j][l] = k;
          }
          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_tran[i][j] = end_wk;
      }
    } else {
      start_wk = 1;
      qc_prob_total = 0.0;

      for (j=0; j<=93; j++) {
        if (uni_ep[i][j] == 0) {
          continue;
        }
        qc_prob_total += uni_ep[i][j];
        end_wk = int(qc_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }
        for (k=start_wk; k<=end_wk; k++) {
          freq2qc[i][k] = j;
        }
        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_freq[i] = end_wk;
    }
  }

  // simulation
  if ((fp = fopen(transcript.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", transcript.file);
    return FAILED;
  }

  if (sim.hp_del_bias == 1) {
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] = 1;
    }
  } else {
    for(i=0; i<=10; i++) {
      transcript.hpfreq[i] = 0;
    }
    flg1 = 1;
    while (fgets(line, BUF_SIZE, fp) != NULL) {
      if (trim(line) == EXISTS_LINE_FEED) {
        flg2 = 1;
      } else {
        flg2 = 0;
      }
      if (flg1 == 1) {
        tp = strtok(line, "\t");
        strncpy(transcript.id, tp, TRANS_ID_LEN_MAX);
        transcript.id[TRANS_ID_LEN_MAX] = '\0';
        tp = strtok(NULL, "\t");
        transcript.plus_exp = atoi(tp);
        tp = strtok(NULL, "\t");
        transcript.minus_exp = atoi(tp);
        tp = strtok(NULL, "\t");
        copy_size = strlen(tp);
        memcpy(transcript.seq, tp, copy_size);
        offset = copy_size;
      } else {
        copy_size = strlen(line);
        memcpy(transcript.seq + offset, line, copy_size);
        offset += copy_size;
      }
      if (flg2 == 1) {
        transcript.seq[offset] = '\0';
        transcript.len = strlen(transcript.seq);
        read_num = transcript.plus_exp + transcript.minus_exp;

        for(i=0; i<transcript.len; i++) {
          transcript.seq[i] = toupper(transcript.seq[i]);
        }

        nstart = 0;
        nend = 0;
        nnum = 1; 
        for(i=1; i<=transcript.len; i++) {
          if ((i < transcript.len) && (transcript.seq[i-1] == transcript.seq[i])) {
            nend = i;
            nnum ++;
            if (nnum > 11) {
              nnum = 10;
            }
          } else {
            if (transcript.seq[i-1] == 'N') {
              transcript.hpfreq[1] += read_num * (nend - nstart + 1);
            } else {
              transcript.hpfreq[nnum] += read_num * (nend - nstart + 1);
            }
            nstart = i;
            nend = nstart;
            nnum = 1;
          }
        }
      }
      flg1 = flg2;
    }
    rewind(fp);

    sum1 = 0;
    sum2 = 0;
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] = 1 + (sim.hp_del_bias - 1) / 9 * (i - 1);
      sum1 += transcript.hpfreq[i] * transcript.hp_del_bias[i];
      sum2 += transcript.hpfreq[i];
    }
    rate = (double)sum2 / sum1;
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] *= rate;
    }
  }

  flg1 = 1;
  while (fgets(line, BUF_SIZE, fp) != NULL) {
    if (trim(line) == EXISTS_LINE_FEED) {
      flg2 = 1;
    } else {
      flg2 = 0;
    }
    if (flg1 == 1) {
      tp = strtok(line, "\t");
      strncpy(transcript.id, tp, TRANS_ID_LEN_MAX);
      transcript.id[TRANS_ID_LEN_MAX] = '\0';
      tp = strtok(NULL, "\t");
      transcript.plus_exp = atoi(tp);
      tp = strtok(NULL, "\t");
      transcript.minus_exp = atoi(tp);
      tp = strtok(NULL, "\t");
      copy_size = strlen(tp);
      memcpy(transcript.seq, tp, copy_size);
      offset = copy_size;
    } else {
      copy_size = strlen(line);
      memcpy(transcript.seq + offset, line, copy_size);
      offset += copy_size;
    }
    if (flg2 == 1) {
      transcript.seq[offset] = '\0';
      transcript.len = strlen(transcript.seq);
      read_num = transcript.plus_exp + transcript.minus_exp;

      for(i=0; i<transcript.len; i++) {
        transcript.seq[i] = toupper(transcript.seq[i]);
      }

      nstart = 0;
      nend = 0;
      nnum = 1; 
      for(i=1; i<=transcript.len; i++) {
        if ((i < transcript.len) && (transcript.seq[i-1] == transcript.seq[i])) {
          nend = i;
          nnum ++;
          if (nnum > 11) {
            nnum = 10;
          }
        } else {
          if (transcript.seq[i-1] == 'N') {
            for(j=nstart; j<=nend; j++) {
              transcript.hp[j] = 1;
            }
          } else {
            for(j=nstart; j<=nend; j++) {
              transcript.hp[j] = nnum;
            }
          }
          nstart = i;
          nend = nstart;
          nnum = 1;
        }
      }

      for (i=1; i<=read_num; i++) {
        index = rand() % len_rand_value + 1;
        mut.len = prob2len[index];

        index = rand() % accuracy_rand_value + 1;
        mut.acc = prob2accuracy[index];

        rank = ceil((double)transcript.len / 1000);
        index = rand() % ssp_rand_value[rank] + 1;
        if (prob2ssp[rank][index] == 0) {
          value = 0.0;
        } else {
          value = ((double)prob2ssp[rank][index] - 2.5) / 100;
        }
        mut.offset = int((double)transcript.len * value + 0.5);
        if (mut.offset + mut.len > transcript.len) {
          mut.len = transcript.len - mut.offset;
        }

        mut.seq_left = mut.offset + 1;
        mut.seq_right = mut.offset + mut.len;
        sim.res_num ++;
        //fprintf(stderr,"%d : %d : %d : %d : %d : %d : %d\n",sim.res_num,transcript.len,rank,mut.offset,mut.len,mut.seq_left,mut.seq_right);

        for (j=0; j<mut.len; j++) {
          mut.seq[j] = transcript.seq[mut.offset + j];
          mut.hp[j] = transcript.hp[mut.offset + j];
        }
        mut.seq[mut.len] = '\0';
        if (i <= transcript.plus_exp) {
          mut.seq_strand = '+';
        } else {
          mut.seq_strand = '-';
          revcomp(mut.seq);
          revshort(mut.hp, mut.len);
        }

        for (h=0; h<sim.pass_num; h++) {
          ref_offset = 0;
          read_offset = 0;
          maf_offset = 0;
          while (ref_offset < mut.len) {
            if (qshmm.exist_hmm[mut.acc] == 1) {
              if (read_offset == 0) {
                index = rand() % qc_rand_value_init[mut.acc] + 1;
                state = init2state[mut.acc][index];
              } else {
                index = rand() % qc_rand_value_tran[mut.acc][state] + 1;
                state = tran2state[mut.acc][state][index];
              }
              index = rand() % qc_rand_value_emis[mut.acc][state] + 1;
              index = emis2qc[mut.acc][state][index];
            } else {
              index = rand() % qc_rand_value_freq[mut.acc] + 1;
              index = freq2qc[mut.acc][index];
            }
            mut.qc[read_offset] = qc[index].character;

            nt = mut.seq[ref_offset];
            qc_value = (int)mut.qc[read_offset] - 33;
            rand_value = rand() % 1000000;
            if (rand_value < mut.sub_thre[qc_value]) {
              sim.res_sub_num ++;
              index = rand() % 3;
              if (nt == 'A') {
                mut.read_seq[read_offset] = mut.sub_nt_a[index];
              } else if (nt == 'T') {
                mut.read_seq[read_offset] = mut.sub_nt_t[index];
              } else if (nt == 'G') {
                mut.read_seq[read_offset] = mut.sub_nt_g[index];
              } else if (nt == 'C') {
                mut.read_seq[read_offset] = mut.sub_nt_c[index];
              } else {
                index = rand() % 4;
                mut.read_seq[read_offset] = mut.sub_nt_n[index];
              }
              mut.maf_ref_seq[maf_offset] = nt;
              ref_offset ++;
            } else if (rand_value < mut.ins_thre[qc_value]) {
              sim.res_ins_num ++;
              index = rand() % 8;
              if (index >= 4) {
                mut.read_seq[read_offset] = nt;
              } else {
                mut.read_seq[read_offset] = mut.ins_nt[index];
              }
              mut.maf_ref_seq[maf_offset] = '-';
            } else {
              mut.read_seq[read_offset] = nt;
              mut.maf_ref_seq[maf_offset] = nt;
              ref_offset ++;
            }
            mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
            maf_offset ++;
            read_offset ++;

            while (ref_offset < mut.len) {
              hp = mut.hp[ref_offset-1];
              rand_value = rand() % 1000000;
              qc_value = (int)mut.qc[read_offset-1] - 33;
              if (rand_value < mut.del_thre[qc_value] * transcript.hp_del_bias[hp]) {
                sim.res_del_num ++;
                mut.maf_seq[maf_offset] = '-';
                mut.maf_ref_seq[maf_offset] = mut.seq[ref_offset];
                maf_offset ++;
                ref_offset ++;
              } else {
                break;
              }
            }
          }
          mut.qc[read_offset] = '\0';
          mut.read_seq[read_offset] = '\0';
          mut.maf_seq[maf_offset] = '\0';
          mut.maf_ref_seq[maf_offset] = '\0';

          if (mut.seq_strand == '-') {
            revcomp(mut.maf_seq);
            revcomp(mut.maf_ref_seq);
          }

          len = strlen(mut.read_seq);
          sim.res_len_total += len;

          freq_len[len] ++;

          if (len > sim.res_len_max) {
            sim.res_len_max = len;
          }
          if (len < sim.res_len_min) {
            sim.res_len_min = len;
          }

          prob = 0.0;
          for (j=0; j<len; j++) {
            prob += qc[(int)mut.qc[j] - 33].prob;
          }
          value = 1.0 - (prob / len);
          accuracy_total += value;
          acc_wk = (int)(value * 100000 + 0.5);
          freq_accuracy[acc_wk] ++;

          if (sim.pass_num == 1) {
            sprintf(id, "%s_%ld", sim.id_prefix, sim.res_num);
            fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.qc);
          } else {
            sprintf(id, "%s/%ld/%ld", sim.id_prefix, sim.res_num, h);
            fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.qc);
            fprintf(fp_sam, "\tcx:i:3\tip:B:C");
            for (j=0; j<len; j++) {
              fprintf(fp_sam, ",9");
            }
            fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
            for (j=0; j<len; j++) {
              fprintf(fp_sam, ",9");
            }
            qeval = len - 1;
            fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
          }

          digit_num1[0] = strlen(transcript.id);
          digit_num2[0] = 1 + count_digit(sim.res_num);
          digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

          digit_num1[1] = count_digit((mut.seq_left - 1));
          digit_num2[1] = 1;
          digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

          digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
          digit_num2[2] = count_digit(len);
          digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

          digit_num1[3] = count_digit(transcript.len);
          digit_num2[3] = count_digit(len);
          digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

          fprintf(fp_maf, "a\ns %s", transcript.id);
          while (digit_num1[0] ++ < digit_num[0]) {
            fprintf(fp_maf, " ");
          }
          while (digit_num1[1] ++ < digit_num[1]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld", mut.seq_left - 1);
          while (digit_num1[2] ++ < digit_num[2]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
          while (digit_num1[3] ++ < digit_num[3]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %s\n", transcript.len, mut.maf_ref_seq);
          fprintf(fp_maf, "s %s", id);
          while (digit_num2[0] ++ < digit_num[0]) {
            fprintf(fp_maf, " ");
          }
          while (digit_num2[1] ++ < digit_num[1]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %d", 0);
          while (digit_num2[2] ++ < digit_num[2]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
          while (digit_num2[3] ++ < digit_num[3]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
        }
      }
    }
    flg1 = flg2;
  }
  fclose(fp);

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}

//////////////////////////////////////////////////////////////
// Function: simulate_by_qshmm_templ - Simulate by Template //
//////////////////////////////////////////////////////////////

int simulate_by_qshmm_templ() {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;
  char *ret_pointer;
  long len;
  long h, i, j, k, l;
  double prob, mean, variance, sd;
  double freq_total, accuracy_prob_total, qc_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;
  static long prob2accuracy[100001];
  static long freq2qc[ACCURACY_MAX+1][1001];
  long state;
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2qc[ACCURACY_MAX+1][STATE_MAX+1][101];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][101];
  long accuracy_rand_value;
  long qc_rand_value_freq[ACCURACY_MAX+1];
  long qc_rand_value_init[ACCURACY_MAX+1];
  long qc_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long qc_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, pre_index;
  long accuracy_min, accuracy_max;
  char id[128], nt;
  int digit_num1[4], digit_num2[4], digit_num[4];
  int qeval;
  long read_offset, ref_offset, maf_offset;
  long qc_value, rand_value;
  int hp;
  long nstart, nend;
  short nnum;
  long sum1, sum2;
  double rate;

  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // quality code distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (qshmm.exist_hmm[i] == 1) {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=1; j<=STATE_MAX; j++) {
        if (qshmm.ip[i][j] == 0) {
          continue;
        }
        qc_prob_total += qshmm.ip[i][j];
        end_wk = int(qc_prob_total * 100 + 0.5);
        if (end_wk > 100) {
          end_wk = 100;
        }

        for (k=start_wk; k<=end_wk; k++) {
          init2state[i][k] = j;
        }

        if (end_wk >= 100) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_init[i] = end_wk;

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=0; k<=93; k++) {
          if (qshmm.ep[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.ep[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            emis2qc[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_emis[i][j] = end_wk;
      }

      for (j=1; j<=STATE_MAX; j++) {
        start_wk = 1; 
        qc_prob_total = 0.0;

        for (k=1; k<=STATE_MAX; k++) {
          if (qshmm.tp[i][j][k] == 0) {
            continue;
          }
          qc_prob_total += qshmm.tp[i][j][k];
          end_wk = int(qc_prob_total * 100 + 0.5);
          if (end_wk > 100) {
            end_wk = 100;
          }

          for (l=start_wk; l<=end_wk; l++) {
            tran2state[i][j][l] = k;
          }

          if (end_wk >= 100) {
            break;
          }
          start_wk = end_wk + 1;
        }
        qc_rand_value_tran[i][j] = end_wk;
      }
    } else {
      start_wk = 1; 
      qc_prob_total = 0.0;

      for (j=0; j<=93; j++) {
        if (uni_ep[i][j] == 0) {
          continue;
        }
        qc_prob_total += uni_ep[i][j];
        end_wk = int(qc_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (k=start_wk; k<=end_wk; k++) {
          freq2qc[i][k] = j;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      qc_rand_value_freq[i] = end_wk;
    }
  }

  // simulation
  if ((fp = fopen(templ.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", templ.file);
    return FAILED;
  }

  if (sim.hp_del_bias == 1) {
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] = 1;
    }
  } else {
    for(i=0; i<=10; i++) {
      templ.hpfreq[i] = 0;
    }
    while (1) {
      ret_pointer = fgets(line, BUF_SIZE, fp);
      if (((ret_pointer == NULL) || (line[0] == '>')) && (offset != 0)) {
        templ.seq[offset] = '\0';
        templ.len = strlen(templ.seq);

        for(i=1; i<=templ.len; i++) {
          templ.seq[i] = toupper(templ.seq[i]);
        }

        nstart = 0;
        nend = 0;
        nnum = 1;
        for(i=1; i<=templ.len; i++) {
          if ((i < templ.len) && (templ.seq[i-1] == templ.seq[i])) {
            nend = i;
            nnum ++;
            if (nnum > 11) {
              nnum = 10;
            }
          } else {
            if (templ.seq[i-1] == 'N') {
              templ.hpfreq[1] += nend - nstart + 1;
            } else {
              templ.hpfreq[nnum] += nend - nstart + 1;
            }
            nstart = i;
            nend = nstart;
            nnum = 1;
          }
        }
      }

      if (ret_pointer == NULL) {
        break;
      }

      ret = trim(line);

      if (line[0] == '>') {
        strncpy(templ.id, line+1, REF_ID_LEN_MAX);
        templ.id[REF_ID_LEN_MAX] = '\0';
        offset = 0;
        while (ret != EXISTS_LINE_FEED) {
          if (fgets(line, BUF_SIZE, fp) == NULL) {
            break;
          }
          ret = trim(line);
        }
      } else {
        copy_size = strlen(line);
        memcpy(templ.seq + offset, line, copy_size);
        offset += copy_size;
      }
    }
    rewind(fp);

    sum1 = 0;
    sum2 = 0;
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] = 1 + (sim.hp_del_bias - 1) / 9 * (i - 1);
      sum1 += templ.hpfreq[i] * templ.hp_del_bias[i];
      sum2 += templ.hpfreq[i];
    }
    rate = (double)sum2 / sum1;
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] *= rate;
    }
  }

  offset = 0;
  while (1) {
    ret_pointer = fgets(line, BUF_SIZE, fp);
    if (((ret_pointer == NULL) || (line[0] == '>')) && (offset != 0)) {
      templ.seq[offset] = '\0';
      templ.len = strlen(templ.seq);

      for(i=1; i<=templ.len; i++) {
        templ.seq[i] = toupper(templ.seq[i]);
      }

      nstart = 0;
      nend = 0;
      nnum = 1;
      for(i=1; i<=templ.len; i++) {
        if ((i < templ.len) && (templ.seq[i-1] == templ.seq[i])) {
          nend = i;
          nnum ++;
          if (nnum > 11) {
            nnum = 10;
          }
        } else {
          if (templ.seq[i-1] == 'N') {
            for(j=nstart; j<=nend; j++) {
              templ.hp[j] = 1;
            }
          } else {
            for(j=nstart; j<=nend; j++) {
              templ.hp[j] = nnum;
            }
          }
          nstart = i;
          nend = nstart;
          nnum = 1;
        }
      }

      index = rand() % accuracy_rand_value + 1;
      mut.acc = prob2accuracy[index];
      mut.seq_left = 1;
      mut.seq_right = templ.len;
      sim.res_num ++;
      mut.seq_strand = '+';
      //fprintf(stderr,"%d : %d : %d : %d\n",sim.res_num,templ.len,mut.seq_left,mut.seq_right);

      for (h=0; h<sim.pass_num; h++) {
        ref_offset = 0;
        read_offset = 0;
        maf_offset = 0;
        while (ref_offset < templ.len) {
          if (qshmm.exist_hmm[mut.acc] == 1) {
            if (read_offset == 0) {
              index = rand() % qc_rand_value_init[mut.acc] + 1;
              state = init2state[mut.acc][index];
            } else {
              index = rand() % qc_rand_value_tran[mut.acc][state] + 1;
              state = tran2state[mut.acc][state][index];
            }
            index = rand() % qc_rand_value_emis[mut.acc][state] + 1;
            index = emis2qc[mut.acc][state][index];
          } else {
            index = rand() % qc_rand_value_freq[mut.acc] + 1;
            index = freq2qc[mut.acc][index];
          }
          mut.qc[read_offset] = qc[index].character;

          nt = templ.seq[ref_offset];
          qc_value = (int)mut.qc[read_offset] - 33;
          rand_value = rand() % 1000000;
          if (rand_value < mut.sub_thre[qc_value]) {
            sim.res_sub_num ++;
            index = rand() % 3;
            if (nt == 'A') {
              mut.read_seq[read_offset] = mut.sub_nt_a[index];
            } else if (nt == 'T') {
              mut.read_seq[read_offset] = mut.sub_nt_t[index];
            } else if (nt == 'G') {
              mut.read_seq[read_offset] = mut.sub_nt_g[index];
            } else if (nt == 'C') {
              mut.read_seq[read_offset] = mut.sub_nt_c[index];
            } else {
              index = rand() % 4;
              mut.read_seq[read_offset] = mut.sub_nt_n[index];
            }
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
          } else if (rand_value < mut.ins_thre[qc_value]) {
            sim.res_ins_num ++;
            index = rand() % 8;
            if (index >= 4) {
              mut.read_seq[read_offset] = nt;
            } else {
              mut.read_seq[read_offset] = mut.ins_nt[index];
            }
            mut.maf_ref_seq[maf_offset] = '-';
          } else {
            mut.read_seq[read_offset] = nt;
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
          }
          mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
          maf_offset ++;
          read_offset ++;

          while (ref_offset < templ.len) {
            hp = templ.hp[ref_offset-1];
            rand_value = rand() % 1000000;
            qc_value = (int)mut.qc[read_offset-1] - 33;
            if (rand_value < mut.del_thre[qc_value] * templ.hp_del_bias[hp]) {
              sim.res_del_num ++;
              mut.maf_seq[maf_offset] = '-';
              mut.maf_ref_seq[maf_offset] = templ.seq[ref_offset];
              maf_offset ++;
              ref_offset ++;
            } else {
              break;
            }
          }
        }
        mut.qc[read_offset] = '\0';
        mut.read_seq[read_offset] = '\0';
        mut.maf_seq[maf_offset] = '\0';
        mut.maf_ref_seq[maf_offset] = '\0';

        len = strlen(mut.read_seq);
        sim.res_len_total += len;
        freq_len[len] ++;

        if (len > sim.res_len_max) {
          sim.res_len_max = len;
        }
        if (len < sim.res_len_min) {
          sim.res_len_min = len;
        }

        prob = 0.0;
        for (i=0; i<len; i++) {
          prob += qc[(int)mut.qc[i] - 33].prob;
        }
        value = 1.0 - (prob / len);
        accuracy_total += value;
        accuracy = (int)(value * 100000 + 0.5);
        freq_accuracy[accuracy] ++;

        if (sim.pass_num == 1) {
          sprintf(id, "%s_%ld", sim.id_prefix, sim.res_num);
          fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.qc);
        } else {
          sprintf(id, "%s/%ld/%ld", sim.id_prefix, sim.res_num, h);
          fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.qc);
          fprintf(fp_sam, "\tcx:i:3\tip:B:C");
          for (i=0; i<len; i++) {
            fprintf(fp_sam, ",9");
          }
          fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
          for (i=0; i<len; i++) {
            fprintf(fp_sam, ",9");
          }
          qeval = len - 1;
          fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
        }

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim.res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut.seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(templ.len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns %s", templ.id);
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut.seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", templ.len, mut.maf_ref_seq);
        fprintf(fp_maf, "s %s", id);
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
      }
    }

    if (ret_pointer == NULL) {
      break;
    }

    ret = trim(line);

    if (line[0] == '>') {
      strncpy(templ.id, line+1, REF_ID_LEN_MAX);
      templ.id[REF_ID_LEN_MAX] = '\0';
      offset = 0;
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(templ.seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}


//////////////////////////////////////////////////////
// Function: simulate_by_errhmm - Simulate by Model //
//////////////////////////////////////////////////////

int simulate_by_errhmm() {
  long len;
  long long len_total = 0;
  long i, j, k, l, h;
  long state;
  double prob, mean, variance, sd;
  double kappa, theta, gamma;
  double len_prob_total, freq_total, accuracy_prob_total, err_prob_total;
  double accuracy_total = 0.0;
  double value, sum;
  static long prob2len[100001], prob2accuracy[100001];
  static long freq2err[ACCURACY_MAX+1][1001];
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2err[ACCURACY_MAX+1][STATE_MAX+1][1001];
  static long emis2del[ACCURACY_MAX+1][STATE_MAX+1];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][1001];
  long len_rand_value, accuracy_rand_value;
  long err_rand_value_freq[ACCURACY_MAX+1];
  long err_rand_value_init[ACCURACY_MAX+1];
  long err_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long err_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, index2, pre_index;
  char id[128];
  char nt, tmperr;
  int acc_wk, accuracy_min, accuracy_max;
  int digit_num1[4], digit_num2[4], digit_num[4];
  int acc_tmp, rate_mag;
  int qeval;
  long read_offset, ref_offset, maf_offset;
  int hp;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // length distribution
  variance = pow(sim.len_sd, 2);
  kappa = pow(sim.len_mean, 2) / variance; 
  theta = variance / sim.len_mean;
  gamma = tgamma(kappa);

  if (sim.len_sd == 0.0) {
    prob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1; 
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += pow(i, kappa-1) * exp(-1 * i / theta) / pow(theta, kappa) / gamma;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }

      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }

      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (sim.pass_num == 1) {
    if (len_rand_value < 1) {
      fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
      return FAILED;
    }
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // error distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (errhmm.exist_hmm[i] == 0) {
      continue;
    }

    start_wk = 1; 
    err_prob_total = 0.0;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      if (errhmm.ip[i][j] == 0) {
        continue;
      }
      err_prob_total += errhmm.ip[i][j];
      end_wk = int(err_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }
      for (k=start_wk; k<=end_wk; k++) {
        init2state[i][k] = j;
      }
      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    err_rand_value_init[i] = end_wk;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1; 
      err_prob_total = 0.0;

      emis2del[i][j] = int(errhmm.ep[i][j][3] * 1000 + 0.5);
      for (k=0; k<=2; k++) {
        if (errhmm.ep[i][j][k] <= 0) {
          continue;
        }
        err_prob_total += errhmm.ep[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (l=start_wk; l<=end_wk; l++) {
          emis2err[i][j][l] = k;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_emis[i][j] = end_wk;
    }

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1; 
      err_prob_total = 0.0;

      for (k=1; k<=STATE_MAX; k++) {
        if (errhmm.tp[i][j][k] == 0) {
          continue;
        }
        err_prob_total += errhmm.tp[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (l=start_wk; l<=end_wk; l++) {
          tran2state[i][j][l] = k;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_tran[i][j] = end_wk;
    }
  }

  // simulation
  while (len_total < sim.len_quota) {
    index = rand() % len_rand_value + 1;
    mut.len = prob2len[index];
    if (len_total + mut.len > sim.len_quota) {
      mut.len = sim.len_quota - len_total;
      if (mut.len < sim.len_min) {
        mut.len = sim.len_min;
      }
    }

    index = rand() % accuracy_rand_value + 1;
    mut.acc = prob2accuracy[index];
    if (mut.len >= genome.len) {
      mut.offset = 0;
      mut.len = genome.len;
    } else {
      mut.offset = rand() % (genome.len - mut.len + 1);
    }

    mut.seq_left = mut.offset + 1;
    mut.seq_right = mut.offset + mut.len;
    sim.res_num ++;

    for (i=0; i<mut.len; i++) {
      mut.seq[i] = genome.seq[mut.offset + i];
      mut.hp[i] = genome.hp[mut.offset + i];
    }
    mut.seq[mut.len] = '\0';
    if (sim.res_num % 2 == 1) {
      mut.seq_strand = '+';
    } else {
      mut.seq_strand = '-';
      revcomp(mut.seq);
      revshort(mut.hp, mut.len);
    }

    //mut.acc = 95;
    if (mut.acc < errhmm.acc_min) {
      rate_mag = int((double)(errhmm.acc_min-mut.acc)/errhmm.acc_min*100);
    } else if (mut.acc > errhmm.acc_max) {
      rate_mag = int((double)(mut.acc-errhmm.acc_max)/(100-errhmm.acc_max)*100);
    }

    for (h=0; h<sim.pass_num; h++) {
      mut.err_num = 0;
      if (mut.acc == 100) {
        for (i=0; i<mut.len; i++) {
          nt = mut.seq[i];
          mut.read_seq[i] = nt;
          mut.maf_ref_seq[i] = nt;
          mut.maf_seq[i] = nt;
        }
        read_offset = mut.len;
        maf_offset = mut.len;
      } else {
        ref_offset = 0;
        read_offset = 0;
        maf_offset = 0;
        while (ref_offset < mut.len) {
          nt = mut.seq[ref_offset];
          if (errhmm.exist_hmm[mut.acc] == 1) {
            if (read_offset == 0) {
              index = rand() % err_rand_value_init[mut.acc] + 1;
              state = init2state[mut.acc][index];
            } else {
              index = rand() % err_rand_value_tran[mut.acc][state] + 1;
              state = tran2state[mut.acc][state][index];
            }
            hp = mut.hp[ref_offset];
            index = rand() % 1000 + 1;
            if (index <= emis2del[mut.acc][state] * genome.hp_del_bias[hp]) {
              index = 3;
            } else {
              if (err_rand_value_emis[mut.acc][state] == 0) {
                index = rand() % 3;
              } else {
                index = rand() % err_rand_value_emis[mut.acc][state] + 1;
                index = emis2err[mut.acc][state][index];
              }
            }
          } else if (mut.acc < errhmm.acc_min) {
            if (read_offset == 0) {
              index = rand() % err_rand_value_init[errhmm.acc_min] + 1;
              state = init2state[errhmm.acc_min][index];
            } else {
              index = rand() % err_rand_value_tran[errhmm.acc_min][state] + 1;
              state = tran2state[errhmm.acc_min][state][index];
            }
            hp = mut.hp[ref_offset];
            index = rand() % 1000 + 1;
            if (index <= emis2del[errhmm.acc_min][state] * genome.hp_del_bias[hp]) {
              index = 3;
            } else {
              if (err_rand_value_emis[errhmm.acc_min][state] == 0) {
                index = rand() % 3;
              } else {
                index = rand() % err_rand_value_emis[errhmm.acc_min][state] + 1;
                index = emis2err[errhmm.acc_min][state][index];
              }
            }
            if (index == 0) {
              index = rand() % 100 + 1;
              if (index <= rate_mag) {
                index = rand() % 3 + 1;
              } else {
                index = 0;
              }
            }
          } else {
            if (read_offset == 0) {
              index = rand() % err_rand_value_init[errhmm.acc_max] + 1;
              state = init2state[errhmm.acc_max][index];
            } else {
              index = rand() % err_rand_value_tran[errhmm.acc_max][state] + 1;
              state = tran2state[errhmm.acc_max][state][index];
            }
            hp = mut.hp[ref_offset];
            index = rand() % 1000 + 1;
            if (index <= emis2del[errhmm.acc_max][state] * genome.hp_del_bias[hp]) {
              index = 3;
            } else {
              if (err_rand_value_emis[errhmm.acc_max][state] == 0) {
                index = rand() % 3;
              } else {
                index = rand() % err_rand_value_emis[errhmm.acc_max][state] + 1;
                index = emis2err[errhmm.acc_max][state][index];
              }
            }
            if (index != 0) {
              index2 = rand() % 100 + 1;
              if (index2 <= rate_mag) {
                index = 0;
              }
            }
          }
          tmperr = err[index].character;

          if (tmperr == '0') {
            mut.read_seq[read_offset] = nt;
            mut.maf_seq[maf_offset] = nt;
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
            read_offset ++;
          } else if (tmperr == '1') {
            mut.err_num ++;
            sim.res_sub_num ++;
            index = rand() % 3;
            if (nt == 'A') {
              mut.read_seq[read_offset] = mut.sub_nt_a[index];
            } else if (nt == 'T') {
              mut.read_seq[read_offset] = mut.sub_nt_t[index];
            } else if (nt == 'G') {
              mut.read_seq[read_offset] = mut.sub_nt_g[index];
            } else if (nt == 'C') {
              mut.read_seq[read_offset] = mut.sub_nt_c[index];
            } else {
              index = rand() % 4;
              mut.read_seq[read_offset] = mut.sub_nt_n[index];
            }
            mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
            read_offset ++;
          } else if (tmperr == '2') {
            mut.err_num ++;
            sim.res_ins_num ++;
            index = rand() % 8;
            if (index >= 4) {
              mut.read_seq[read_offset] = nt;
            } else {
              mut.read_seq[read_offset] = mut.ins_nt[index];
            }
            mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
            mut.maf_ref_seq[maf_offset] = '-';
            read_offset ++;
          } else {
            mut.err_num ++;
            sim.res_del_num ++;
            mut.maf_seq[maf_offset] = '-';
            mut.maf_ref_seq[maf_offset] = nt;
            ref_offset ++;
          }
          maf_offset ++;
        }
      }
      mut.read_seq[read_offset] = '\0';
      mut.maf_seq[maf_offset] = '\0';
      mut.maf_ref_seq[maf_offset] = '\0';

      if (mut.seq_strand == '-') {
        revcomp(mut.maf_seq);
        revcomp(mut.maf_ref_seq);
      }

      len = strlen(mut.read_seq);
      sim.res_len_total += len;

      if (h == 0) {
        len_total += len;
      }

      freq_len[len] ++;

      if (len > sim.res_len_max) {
        sim.res_len_max = len;
      }
      if (len < sim.res_len_min) {
        sim.res_len_min = len;
      }
   
      value = 1.0 - ((double)mut.err_num / len);
      accuracy_total += value;
      acc_wk = (int)(value * 100000 + 0.5);
      freq_accuracy[acc_wk] ++;

      for (i=0; i<len; i++) {
        mut.new_qc[i] = '!';
      }
      mut.new_qc[len] = '\0';

      if (sim.pass_num == 1) {
        sprintf(id, "%s%ld_%ld", sim.id_prefix, genome.num, sim.res_num);
        fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.new_qc);
      } else {
        sprintf(id, "%s%ld/%ld/%ld", sim.id_prefix, genome.num, sim.res_num, h);
        fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.new_qc);
        fprintf(fp_sam, "\tcx:i:3\tip:B:C");
        for (i=0; i<len; i++) {
          fprintf(fp_sam, ",9");
        }
        fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
        for (i=0; i<len; i++) {
          fprintf(fp_sam, ",9");
        }
        qeval = len - 1;
        fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
      }

      digit_num1[0] = 3;
      digit_num2[0] = 1 + count_digit(sim.res_num);
      digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

      digit_num1[1] = count_digit((mut.seq_left - 1));
      digit_num2[1] = 1;
      digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

      digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
      digit_num2[2] = count_digit(len);
      digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

      digit_num1[3] = count_digit(genome.len);
      digit_num2[3] = count_digit(len);
      digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

      digit_num1[0] = 3;
      fprintf(fp_maf, "a\ns ref");
      while (digit_num1[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num1[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld", mut.seq_left - 1);
      while (digit_num1[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
      while (digit_num1[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n", genome.len, mut.maf_ref_seq);
      fprintf(fp_maf, "s %s", id);
      while (digit_num2[0] ++ < digit_num[0]) {
        fprintf(fp_maf, " ");
      }
      while (digit_num2[1] ++ < digit_num[1]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %d", 0);
      while (digit_num2[2] ++ < digit_num[2]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
      while (digit_num2[3] ++ < digit_num[3]) {
        fprintf(fp_maf, " ");
      }
      fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
    }
  }

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}

////////////////////////////////////////////////////////////
// Function: simulate_by_errhmm_trans - Simulate by Model //
////////////////////////////////////////////////////////////

int simulate_by_errhmm_trans() {
  FILE *fp;
  char *tp;
  char line[BUF_SIZE];
  int flg1, flg2;
  long offset, copy_size;
  long len;
  long i, j, k, l, h;
  long state;
  double prob, mean, variance, sd;
  double kappa, theta, gamma;
  double len_prob_total, freq_total, accuracy_prob_total, err_prob_total;
  double ssp_prob_total;
  double value, sum;
  double accuracy_total = 0.0;
  static long prob2len[100001], prob2accuracy[100001];
  static long prob2ssp[TR_RANK_MAX][1001];
  static long freq2err[ACCURACY_MAX+1][1001];
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2err[ACCURACY_MAX+1][STATE_MAX+1][1001];
  static long emis2del[ACCURACY_MAX+1][STATE_MAX+1];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][1001];
  long len_rand_value, accuracy_rand_value;
  long ssp_rand_value[TR_RANK_MAX+1];
  long err_rand_value_freq[ACCURACY_MAX+1];
  long err_rand_value_init[ACCURACY_MAX+1];
  long err_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long err_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, index2, pre_index;
  char id[128];
  int acc_wk, accuracy_min, accuracy_max;
  int digit_num1[4], digit_num2[4], digit_num[4];
  int acc_tmp, rate_mag;
  int qeval;
  int read_num;
  int rank;
  long read_offset, ref_offset, maf_offset;
  int hp;
  char nt, tmperr;
  long nstart, nend;
  short nnum;
  double rate;
  long sum1, sum2;

  for (i=0; i<=sim.len_max; i++) {
    freq_len[i] = 0;
  }
  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // length distribution
  variance = pow(sim.len_sd, 2);
  kappa = pow(sim.len_mean, 2) / variance; 
  theta = variance / sim.len_mean;
  gamma = tgamma(kappa);

  if (sim.len_sd == 0.0) {
    prob2len[1] = int(sim.len_mean + 0.5);
    len_rand_value = 1;
  } else {
    start_wk = 1; 
    len_prob_total = 0.0;
    for (i=sim.len_min; i<=sim.len_max; i++) {
      len_prob_total += pow(i, kappa-1) * exp(-1 * i / theta) / pow(theta, kappa) / gamma;
      end_wk = int(len_prob_total * 100000 + 0.5);
      if (end_wk > 100000) {
        end_wk = 100000;
      }
      for (j=start_wk; j<=end_wk; j++) {
        prob2len[j] = i;
      }
      if (end_wk >= 100000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    len_rand_value = end_wk;
  }

  if (len_rand_value < 1) {
    fprintf(stderr, "ERROR: length parameters are not appropriate.\n");
    return FAILED;
  }

  // sequencing start pos distribution
  for (i=1; i<=transcript.rank_max; i++) {
    sum = 0;
    value = (double)1 / i;
    for (j=1; j<=21; j++) {
      sum += value / pow(j,(1+value));
    }
    start_wk = 1;
    ssp_prob_total = 0.0;
    for (j=1; j<=21; j++) {
      ssp_prob_total += (value / pow(j,(1+value))) / sum;
      end_wk = int(ssp_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }
      for (k=start_wk; k<=end_wk; k++) {
        prob2ssp[i][k] = (j-1) * 5;
      }
      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    ssp_rand_value[i] = end_wk;
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // error distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (errhmm.exist_hmm[i] == 0) {
      continue;
    }

    start_wk = 1; 
    err_prob_total = 0.0;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      if (errhmm.ip[i][j] == 0) {
        continue;
      }
      err_prob_total += errhmm.ip[i][j];
      end_wk = int(err_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }
      for (k=start_wk; k<=end_wk; k++) {
        init2state[i][k] = j;
      }
      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    err_rand_value_init[i] = end_wk;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1; 
      err_prob_total = 0.0;

      emis2del[i][j] = int(errhmm.ep[i][j][3] * 1000 + 0.5);
      for (k=0; k<=2; k++) {
        if (errhmm.ep[i][j][k] == 0) {
          continue;
        }
        err_prob_total += errhmm.ep[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }
        for (l=start_wk; l<=end_wk; l++) {
          emis2err[i][j][l] = k;
        }
        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_emis[i][j] = end_wk;
    }

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1; 
      err_prob_total = 0.0;

      for (k=1; k<=STATE_MAX; k++) {
        if (errhmm.tp[i][j][k] == 0) {
          continue;
        }
        err_prob_total += errhmm.tp[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }
        for (l=start_wk; l<=end_wk; l++) {
          tran2state[i][j][l] = k;
        }
        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_tran[i][j] = end_wk;
    }
  }

  // simulation
  if ((fp = fopen(transcript.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", transcript.file);
    return FAILED;
  }

  if (sim.hp_del_bias == 1) {
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] = 1;
    }
  } else {
    for(i=0; i<=10; i++) {
      transcript.hpfreq[i] = 0;
    }
    flg1 = 1;
    while (fgets(line, BUF_SIZE, fp) != NULL) {
      if (trim(line) == EXISTS_LINE_FEED) {
        flg2 = 1;
      } else {
        flg2 = 0;
      }
      if (flg1 == 1) {
        tp = strtok(line, "\t");
        strncpy(transcript.id, tp, TRANS_ID_LEN_MAX);
        transcript.id[TRANS_ID_LEN_MAX] = '\0';
        tp = strtok(NULL, "\t");
        transcript.plus_exp = atoi(tp);
        tp = strtok(NULL, "\t");
        transcript.minus_exp = atoi(tp);
        tp = strtok(NULL, "\t");
        copy_size = strlen(tp);
        memcpy(transcript.seq, tp, copy_size);
        offset = copy_size;
      } else {
        copy_size = strlen(line);
        memcpy(transcript.seq + offset, line, copy_size);
        offset += copy_size;
      }
      if (flg2 == 1) {
        transcript.seq[offset] = '\0';
        transcript.len = strlen(transcript.seq);
        read_num = transcript.plus_exp + transcript.minus_exp;

        for(i=1; i<=transcript.len; i++) {
          transcript.seq[i] = toupper(transcript.seq[i]);
        }

        nstart = 0;
        nend = 0;
        nnum = 1;
        for(i=1; i<=transcript.len; i++) {
          if ((i < transcript.len) && (transcript.seq[i-1] == transcript.seq[i])) {
            nend = i;
            nnum ++;
            if (nnum > 11) {
              nnum = 10;
            }
          } else {
            if (transcript.seq[i-1] == 'N') {
              transcript.hpfreq[1] += read_num * (nend - nstart + 1);
            } else {
              transcript.hpfreq[nnum] += read_num * (nend - nstart + 1);
            }
            nstart = i;
            nend = nstart;
            nnum = 1;
          }
        }
      }
      flg1 = flg2;
    }
    rewind(fp);

    sum1 = 0;
    sum2 = 0;
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] = 1 + (sim.hp_del_bias - 1) / 9 * (i - 1);
      sum1 += transcript.hpfreq[i] * transcript.hp_del_bias[i];
      sum2 += transcript.hpfreq[i];
    }
    rate = (double)sum2 / sum1;
    for(i=1; i<=10; i++) {
      transcript.hp_del_bias[i] *= rate;
    }
  }

  flg1 = 1;
  while (fgets(line, BUF_SIZE, fp) != NULL) {
    if (trim(line) == EXISTS_LINE_FEED) {
      flg2 = 1;
    } else {
      flg2 = 0;
    }
    if (flg1 == 1) {
      tp = strtok(line, "\t");
      strncpy(transcript.id, tp, TRANS_ID_LEN_MAX);
      transcript.id[TRANS_ID_LEN_MAX] = '\0';
      tp = strtok(NULL, "\t");
      transcript.plus_exp = atoi(tp);
      tp = strtok(NULL, "\t");
      transcript.minus_exp = atoi(tp);
      tp = strtok(NULL, "\t");
      copy_size = strlen(tp);
      memcpy(transcript.seq, tp, copy_size);
      offset = copy_size;
    } else {
      copy_size = strlen(line);
      memcpy(transcript.seq + offset, line, copy_size);
      offset += copy_size;
    }
    if (flg2 == 1) {
      transcript.seq[offset] = '\0';
      transcript.len = strlen(transcript.seq);
      read_num = transcript.plus_exp + transcript.minus_exp;

      for(i=1; i<=transcript.len; i++) {
        transcript.seq[i] = toupper(transcript.seq[i]);
      }

      nstart = 0;
      nend = 0;
      nnum = 1;
      for(i=1; i<=transcript.len; i++) {
        if ((i < transcript.len) && (transcript.seq[i-1] == transcript.seq[i])) {
          nend = i;
          nnum ++;
          if (nnum > 11) {
            nnum = 10;
          }
        } else {
          if (transcript.seq[i-1] == 'N') {
            for(j=nstart; j<=nend; j++) {
              transcript.hp[j] = 1;
            }
          } else {
            for(j=nstart; j<=nend; j++) {
              transcript.hp[j] = nnum;
            }
          }
          nstart = i;
          nend = nstart;
          nnum = 1;
        }
      }

      for (i=1; i<=read_num; i++) {
        index = rand() % len_rand_value + 1;
        mut.len = prob2len[index];

        index = rand() % accuracy_rand_value + 1;
        mut.acc = prob2accuracy[index];

        rank = ceil((double)transcript.len / 1000);
        index = rand() % ssp_rand_value[rank] + 1;
        if (prob2ssp[rank][index] == 0) {
          value = 0.0;
        } else {
          value = ((double)prob2ssp[rank][index] - 2.5) / 100;
        }
        mut.offset = int((double)transcript.len * value + 0.5);
        if (mut.offset + mut.len > transcript.len) {
          mut.len = transcript.len - mut.offset;
        }

        mut.seq_left = mut.offset + 1;
        mut.seq_right = mut.offset + mut.len;
        sim.res_num ++;
        //fprintf(stderr,"%d : %d : %d : %d : %d : %d : %d\n",sim.res_num,transcript.len,rank,mut.offset,mut.len,mut.seq_left,mut.seq_right);
 
        for (j=0; j<mut.len; j++) {
          mut.seq[j] = transcript.seq[mut.offset + j];
          mut.hp[j] = transcript.hp[mut.offset + j];
        }
        mut.seq[mut.len] = '\0';
        if (i <= transcript.plus_exp) {
          mut.seq_strand = '+';
        } else {
          mut.seq_strand = '-';
          revcomp(mut.seq);
          revshort(mut.hp, mut.len);
        }

        if (mut.acc < errhmm.acc_min) {
          rate_mag = int((double)(errhmm.acc_min-mut.acc)/errhmm.acc_min*100);
        } else if (mut.acc > errhmm.acc_max) {
          rate_mag = int((double)(mut.acc-errhmm.acc_max)/(100-errhmm.acc_max)*100);
        }

        for (h=0; h<sim.pass_num; h++) {
          mut.err_num = 0;
          if (mut.acc == 100) {
            for (i=0; i<mut.len; i++) {
              nt = mut.seq[i];
              mut.read_seq[i] = nt;
              mut.maf_ref_seq[i] = nt;
              mut.maf_seq[i] = nt;
            }
            read_offset = mut.len;
            maf_offset = mut.len;
          } else {
            ref_offset = 0;
            read_offset = 0;
            maf_offset = 0;
            while (ref_offset < mut.len) {
              nt = mut.seq[ref_offset];
              if (errhmm.exist_hmm[mut.acc] == 1) {
                if (read_offset == 0) {
                  index = rand() % err_rand_value_init[mut.acc] + 1;
                  state = init2state[mut.acc][index];
                } else {
                  index = rand() % err_rand_value_tran[mut.acc][state] + 1;
                  state = tran2state[mut.acc][state][index];
                }
                hp = mut.hp[ref_offset];
                index = rand() % 1000 + 1;
                if (index <= emis2del[mut.acc][state] * transcript.hp_del_bias[hp]) {
                  index = 3;
                } else {
                  if (err_rand_value_emis[mut.acc][state] == 0) {
                    index = rand() % 3;
                  } else {
                    index = rand() % err_rand_value_emis[mut.acc][state] + 1;
                    index = emis2err[mut.acc][state][index];
                  }
                }
              } else if (mut.acc < errhmm.acc_min) {
                if (read_offset == 0) {
                  index = rand() % err_rand_value_init[errhmm.acc_min] + 1;
                  state = init2state[errhmm.acc_min][index];
                } else {
                  index = rand() % err_rand_value_tran[errhmm.acc_min][state] + 1;
                  state = tran2state[errhmm.acc_min][state][index];
                }
                hp = mut.hp[ref_offset];
                index = rand() % 1000 + 1;
                if (index <= emis2del[errhmm.acc_min][state] * transcript.hp_del_bias[hp]) {
                  index = 3;
                } else {
                  if (err_rand_value_emis[errhmm.acc_min][state] == 0) {
                    index = rand() % 3;
                  } else {
                    index = rand() % err_rand_value_emis[errhmm.acc_min][state] + 1;
                    index = emis2err[errhmm.acc_min][state][index];
                  }
                }
                if (index == 0) {
                  index = rand() % 100 + 1;
                  if (index <= rate_mag) {
                    index = rand() % 3 + 1;
                  } else {
                    index = 0;
                  }
                }
              } else {
                if (read_offset == 0) {
                  index = rand() % err_rand_value_init[errhmm.acc_max] + 1;
                  state = init2state[errhmm.acc_max][index];
                } else {
                  index = rand() % err_rand_value_tran[errhmm.acc_max][state] + 1;
                  state = tran2state[errhmm.acc_max][state][index];
                }
                hp = mut.hp[ref_offset];
                index = rand() % 1000 + 1;
                if (index <= emis2del[errhmm.acc_max][state] * transcript.hp_del_bias[hp]) {
                  index = 3;
                } else {
                  if (err_rand_value_emis[errhmm.acc_max][state] == 0) {
                    index = rand() % 3;
                  } else {
                    index = rand() % err_rand_value_emis[errhmm.acc_max][state] + 1;
                    index = emis2err[errhmm.acc_max][state][index];
                  }
                }
                if (index != 0) {
                  index2 = rand() % 100 + 1;
                  if (index2 <= rate_mag) {
                    index = 0;
                  }
                }
              }
              tmperr = err[index].character;

              if (tmperr == '0') {
                mut.read_seq[read_offset] = nt;
                mut.maf_seq[maf_offset] = nt;
                mut.maf_ref_seq[maf_offset] = nt;
                ref_offset ++;
                read_offset ++;
              } else if (tmperr == '1') {
                mut.err_num ++;
                sim.res_sub_num ++;
                index = rand() % 3;
                if (nt == 'A') {
                  mut.read_seq[read_offset] = mut.sub_nt_a[index];
                } else if (nt == 'T') {
                  mut.read_seq[read_offset] = mut.sub_nt_t[index];
                } else if (nt == 'G') {
                  mut.read_seq[read_offset] = mut.sub_nt_g[index];
                } else if (nt == 'C') {
                  mut.read_seq[read_offset] = mut.sub_nt_c[index];
                } else {
                  index = rand() % 4;
                  mut.read_seq[read_offset] = mut.sub_nt_n[index];
                }
                mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
                mut.maf_ref_seq[maf_offset] = nt;
                ref_offset ++;
                read_offset ++;
              } else if (tmperr == '2') {
                mut.err_num ++;
                sim.res_ins_num ++;
                index = rand() % 8;
                if (index >= 4) {
                  mut.read_seq[read_offset] = nt;
                } else {
                  mut.read_seq[read_offset] = mut.ins_nt[index];
                }
                mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
                mut.maf_ref_seq[maf_offset] = '-';
                read_offset ++;
              } else {
                mut.err_num ++;
                sim.res_del_num ++;
                mut.maf_seq[maf_offset] = '-';
                mut.maf_ref_seq[maf_offset] = nt;
                ref_offset ++;
              }
              maf_offset ++;
            }
          }
          mut.read_seq[read_offset] = '\0';
          mut.maf_seq[maf_offset] = '\0';
          mut.maf_ref_seq[maf_offset] = '\0';

          if (mut.seq_strand == '-') {
            revcomp(mut.maf_seq);
            revcomp(mut.maf_ref_seq);
          }

          len = strlen(mut.read_seq);
          sim.res_len_total += len;
          freq_len[len] ++;

          if (len > sim.res_len_max) {
            sim.res_len_max = len;
          }
          if (len < sim.res_len_min) {
            sim.res_len_min = len;
          }
   
          value = 1.0 - ((double)mut.err_num / len);
          accuracy_total += value;
          acc_wk = (int)(value * 100000 + 0.5);
          freq_accuracy[acc_wk] ++;

          for (j=0; j<len; j++) {
            mut.new_qc[j] = '!';
          }
          mut.new_qc[len] = '\0';

          if (sim.pass_num == 1) {
            sprintf(id, "%s_%ld", sim.id_prefix, sim.res_num);
            fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.new_qc);
          } else {
            sprintf(id, "%s/%ld/%ld", sim.id_prefix, sim.res_num, h);
            fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.new_qc);
            fprintf(fp_sam, "\tcx:i:3\tip:B:C");
            for (j=0; j<len; j++) {
              fprintf(fp_sam, ",9");
            }
            fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
            for (j=0; j<len; j++) {
              fprintf(fp_sam, ",9");
            }
            qeval = len - 1;
            fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
          }

          digit_num1[0] = strlen(transcript.id);
          digit_num2[0] = 1 + count_digit(sim.res_num);
          digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

          digit_num1[1] = count_digit((mut.seq_left - 1));
          digit_num2[1] = 1;
          digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

          digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
          digit_num2[2] = count_digit(len);
          digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

          digit_num1[3] = count_digit(transcript.len);
          digit_num2[3] = count_digit(len);
          digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

          fprintf(fp_maf, "a\ns %s", transcript.id);
          while (digit_num1[0] ++ < digit_num[0]) {
            fprintf(fp_maf, " ");
          }
          while (digit_num1[1] ++ < digit_num[1]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld", mut.seq_left - 1);
          while (digit_num1[2] ++ < digit_num[2]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
          while (digit_num1[3] ++ < digit_num[3]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %s\n", transcript.len, mut.maf_ref_seq);
          fprintf(fp_maf, "s %s", id);
          while (digit_num2[0] ++ < digit_num[0]) {
            fprintf(fp_maf, " ");
          }
          while (digit_num2[1] ++ < digit_num[1]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %d", 0);
          while (digit_num2[2] ++ < digit_num[2]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
          while (digit_num2[3] ++ < digit_num[3]) {
            fprintf(fp_maf, " ");
          }
          fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
        }
      }
    }
    flg1 = flg2;
  }
  fclose(fp);

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}

///////////////////////////////////////////////////////////////
// Function: simulate_by_errhmm_templ - Simulate by Template //
///////////////////////////////////////////////////////////////

int simulate_by_errhmm_templ() {
  FILE *fp;
  char line[BUF_SIZE];
  long offset = 0;
  long copy_size;
  int ret;
  char *ret_pointer;
  long len;
  long h, i, j, k, l;
  long state;
  double prob, mean, variance, sd;
  double freq_total, accuracy_prob_total, err_prob_total, value, sum;
  double accuracy_total = 0.0;
  int accuracy;
  static long prob2accuracy[100001];
  static long freq2err[ACCURACY_MAX+1][1001];
  static long init2state[ACCURACY_MAX+1][1001];
  static long emis2err[ACCURACY_MAX+1][STATE_MAX+1][1001];
  static long emis2del[ACCURACY_MAX+1][STATE_MAX+1];
  static long tran2state[ACCURACY_MAX+1][STATE_MAX+1][1001];
  long accuracy_rand_value;
  long err_rand_value_freq[ACCURACY_MAX+1];
  long err_rand_value_init[ACCURACY_MAX+1];
  long err_rand_value_emis[ACCURACY_MAX+1][STATE_MAX+1];
  long err_rand_value_tran[ACCURACY_MAX+1][STATE_MAX+1];
  long start_wk, end_wk;
  long index, index2, pre_index;
  long accuracy_min, accuracy_max;
  char id[128];
  int digit_num1[4], digit_num2[4], digit_num[4];
  int acc_tmp, rate_mag;
  int qeval;
  long read_offset, ref_offset, maf_offset;
  int hp;
  char nt, tmperr;
  long nstart, nend;
  short nnum;
  long sum1, sum2;
  double rate;

  for (i=0; i<=100000; i++) {
    freq_accuracy[i] = 0;
  }

  // accuracy distribution
  mean = sim.accuracy_mean * 100;
  accuracy_max = floor(mean * 1.05);
  accuracy_min = floor(mean * 0.75);
  if (accuracy_max > 100) {
    accuracy_max = 100;
  }

  freq_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    freq_total += exp(0.22 * i);
  }
  start_wk = 1; 
  accuracy_prob_total = 0.0;
  for (i=accuracy_min; i<=accuracy_max; i++) {
    accuracy_prob_total += exp(0.22 * i) / freq_total;
    end_wk = int(accuracy_prob_total * 100000 + 0.5);
    if (end_wk > 100000) {
      end_wk = 100000;
    }

    for (j=start_wk; j<=end_wk; j++) {
      prob2accuracy[j] = i;
    }

    if (end_wk >= 100000) {
      break;
    }
    start_wk = end_wk + 1;
  }
  accuracy_rand_value = end_wk;

  if (accuracy_rand_value < 1) {
    fprintf(stderr, "ERROR: accuracy parameters are not appropriate.\n");
    return FAILED;
  }

  // error distribution
  for (i=accuracy_min; i<=accuracy_max; i++) {

    if (errhmm.exist_hmm[i] == 0) {
      continue;
    }

    start_wk = 1;
    err_prob_total = 0.0;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      if (errhmm.ip[i][j] == 0) {
        continue;
      }
      err_prob_total += errhmm.ip[i][j];
      end_wk = int(err_prob_total * 1000 + 0.5);
      if (end_wk > 1000) {
        end_wk = 1000;
      }
      for (k=start_wk; k<=end_wk; k++) {
        init2state[i][k] = j;
      }
      if (end_wk >= 1000) {
        break;
      }
      start_wk = end_wk + 1;
    }
    err_rand_value_init[i] = end_wk;

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1;
      err_prob_total = 0.0;

      emis2del[i][j] = int(errhmm.ep[i][j][3] * 1000 + 0.5);
      for (k=0; k<=2; k++) {
        if (errhmm.ep[i][j][k] == 0) {
          continue;
        }
        err_prob_total += errhmm.ep[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (l=start_wk; l<=end_wk; l++) {
          emis2err[i][j][l] = k;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_emis[i][j] = end_wk;
    }

    for (j=1; j<=errhmm.state_max[i]; j++) {
      start_wk = 1;
      err_prob_total = 0.0;

      for (k=1; k<=STATE_MAX; k++) {
        if (errhmm.tp[i][j][k] == 0) {
          continue;
        }
        err_prob_total += errhmm.tp[i][j][k];
        end_wk = int(err_prob_total * 1000 + 0.5);
        if (end_wk > 1000) {
          end_wk = 1000;
        }

        for (l=start_wk; l<=end_wk; l++) {
          tran2state[i][j][l] = k;
        }

        if (end_wk >= 1000) {
          break;
        }
        start_wk = end_wk + 1;
      }
      err_rand_value_tran[i][j] = end_wk;
    }
  }

  // simulation
  if ((fp = fopen(templ.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", templ.file);
    return FAILED;
  }

  if (sim.hp_del_bias == 1) {
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] = 1;
    }
  } else {
    for(i=0; i<=10; i++) {
      templ.hpfreq[i] = 0;
    }
    while (1) {
      ret_pointer = fgets(line, BUF_SIZE, fp);
      if (((ret_pointer == NULL) || (line[0] == '>')) && (offset != 0)) {
        templ.seq[offset] = '\0';
        templ.len = strlen(templ.seq);

        for(i=1; i<=templ.len; i++) {
          templ.seq[i] = toupper(templ.seq[i]);
        }

        nstart = 0;
        nend = 0;
        nnum = 1;
        for(i=1; i<=templ.len; i++) {
          if ((i < templ.len) && (templ.seq[i-1] == templ.seq[i])) {
            nend = i;
            nnum ++;
            if (nnum > 11) {
              nnum = 10;
            }
          } else {
            if (templ.seq[i-1] == 'N') {
              templ.hpfreq[1] += nend - nstart + 1;
            } else {
              templ.hpfreq[nnum] += nend - nstart + 1;
            }
            nstart = i;
            nend = nstart;
            nnum = 1;
          }
        }
      }

      if (ret_pointer == NULL) {
        break;
      }

      ret = trim(line);

      if (line[0] == '>') {
        strncpy(templ.id, line+1, REF_ID_LEN_MAX);
        templ.id[REF_ID_LEN_MAX] = '\0';
        offset = 0;
        while (ret != EXISTS_LINE_FEED) {
          if (fgets(line, BUF_SIZE, fp) == NULL) {
            break;
          }
          ret = trim(line);
        }
      } else {
        copy_size = strlen(line);
        memcpy(templ.seq + offset, line, copy_size);
        offset += copy_size;
      }
    }
    rewind(fp);

    sum1 = 0;
    sum2 = 0;
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] = 1 + (sim.hp_del_bias - 1) / 9 * (i - 1);
      sum1 += templ.hpfreq[i] * templ.hp_del_bias[i];
      sum2 += templ.hpfreq[i];
    }
    rate = (double)sum2 / sum1;
    for(i=1; i<=10; i++) {
      templ.hp_del_bias[i] *= rate;
    }
  }

  offset = 0;
  while (1) {
    ret_pointer = fgets(line, BUF_SIZE, fp);
    if (((ret_pointer == NULL) || (line[0] == '>')) && (offset != 0)) {
      templ.seq[offset] = '\0';
      templ.len = strlen(templ.seq);

      for(i=1; i<=templ.len; i++) {
        templ.seq[i] = toupper(templ.seq[i]);
      }

      nstart = 0;
      nend = 0;
      nnum = 1;
      for(i=1; i<=templ.len; i++) {
        if ((i < templ.len) && (templ.seq[i-1] == templ.seq[i])) {
          nend = i;
          nnum ++;
          if (nnum > 11) {
            nnum = 10;
          }
        } else {
          if (templ.seq[i-1] == 'N') {
            for(j=nstart; j<=nend; j++) {
              templ.hp[j] = 1;
            }
          } else {
            for(j=nstart; j<=nend; j++) {
              templ.hp[j] = nnum;
            }
          }
          nstart = i;
          nend = nstart;
          nnum = 1;
        }
      }

      index = rand() % accuracy_rand_value + 1;
      mut.acc = prob2accuracy[index];
      mut.seq_left = 1;
      mut.seq_right = templ.len;
      sim.res_num ++;
      mut.seq_strand = '+';

      if (mut.acc < errhmm.acc_min) {
        rate_mag = int((double)(errhmm.acc_min-mut.acc)/errhmm.acc_min*100);
      } else if (mut.acc > errhmm.acc_max) {
        rate_mag = int((double)(mut.acc-errhmm.acc_max)/(100-errhmm.acc_max)*100);
      }

      for (h=0; h<sim.pass_num; h++) {
        mut.err_num = 0;
        if (mut.acc == 100) {
          for (i=0; i<templ.len; i++) {
            nt = templ.seq[i];
            mut.read_seq[i] = nt;
            mut.maf_ref_seq[i] = nt;
            mut.maf_seq[i] = nt;
          }
          read_offset = templ.len;
          maf_offset = templ.len;
        } else {
          ref_offset = 0;
          read_offset = 0;
          maf_offset = 0;
          while (ref_offset < templ.len) {
            nt = templ.seq[ref_offset];
            if (errhmm.exist_hmm[mut.acc] == 1) {
              if (read_offset == 0) {
                index = rand() % err_rand_value_init[mut.acc] + 1;
                state = init2state[mut.acc][index];
              } else {
                index = rand() % err_rand_value_tran[mut.acc][state] + 1;
                state = tran2state[mut.acc][state][index];
              }
              hp = templ.hp[ref_offset];
              index = rand() % 1000 + 1;
              if (index <= emis2del[mut.acc][state] * templ.hp_del_bias[hp]) {
                index = 3;
              } else {
                if (err_rand_value_emis[mut.acc][state] == 0) {
                  index = rand() % 3;
                } else {
                  index = rand() % err_rand_value_emis[mut.acc][state] + 1;
                  index = emis2err[mut.acc][state][index];
                }
              }
            } else if (mut.acc < errhmm.acc_min) {
              if (read_offset == 0) {
                index = rand() % err_rand_value_init[errhmm.acc_min] + 1;
                state = init2state[errhmm.acc_min][index];
              } else {
                index = rand() % err_rand_value_tran[errhmm.acc_min][state] + 1;
                state = tran2state[errhmm.acc_min][state][index];
              }
              hp = templ.hp[ref_offset];
              index = rand() % 1000 + 1;
              if (index <= emis2del[errhmm.acc_min][state] * templ.hp_del_bias[hp]) {
                index = 3;
              } else {
                if (err_rand_value_emis[errhmm.acc_min][state] == 0) {
                  index = rand() % 3;
                } else {
                  index = rand() % err_rand_value_emis[errhmm.acc_min][state] + 1;
                  index = emis2err[errhmm.acc_min][state][index];
                }
              }
              if (index == 0) {
                index = rand() % 100 + 1;
                if (index <= rate_mag) {
                  index = rand() % 3 + 1;
                } else {
                  index = 0;
                }
              }
            } else {
              if (read_offset == 0) {
                index = rand() % err_rand_value_init[errhmm.acc_max] + 1;
                state = init2state[errhmm.acc_max][index];
              } else {
                index = rand() % err_rand_value_tran[errhmm.acc_max][state] + 1;
                state = tran2state[errhmm.acc_max][state][index];
              }
              hp = templ.hp[ref_offset];
              index = rand() % 1000 + 1;
              if (index <= emis2del[errhmm.acc_max][state] * templ.hp_del_bias[hp]) {
                index = 3;
              } else {
                if (err_rand_value_emis[errhmm.acc_max][state] == 0) {
                  index = rand() % 3;
                } else {
                  index = rand() % err_rand_value_emis[errhmm.acc_max][state] + 1;
                  index = emis2err[errhmm.acc_max][state][index];
                }
              }
              if (index != 0) {
                index2 = rand() % 100 + 1;
                if (index2 <= rate_mag) {
                  index = 0;
                }
              }
            }
            tmperr = err[index].character;

            if (tmperr == '0') {
              mut.read_seq[read_offset] = nt;
              mut.maf_seq[maf_offset] = nt;
              mut.maf_ref_seq[maf_offset] = nt;
              ref_offset ++;
              read_offset ++;
            } else if (tmperr == '1') {
              mut.err_num ++;
              sim.res_sub_num ++;
              index = rand() % 3;
              if (nt == 'A') {
                mut.read_seq[read_offset] = mut.sub_nt_a[index];
              } else if (nt == 'T') {
                mut.read_seq[read_offset] = mut.sub_nt_t[index];
              } else if (nt == 'G') {
                mut.read_seq[read_offset] = mut.sub_nt_g[index];
              } else if (nt == 'C') {
                mut.read_seq[read_offset] = mut.sub_nt_c[index];
              } else {
                index = rand() % 4;
                mut.read_seq[read_offset] = mut.sub_nt_n[index];
              }
              mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
              mut.maf_ref_seq[maf_offset] = nt;
              ref_offset ++;
              read_offset ++;
            } else if (tmperr == '2') {
              mut.err_num ++;
              sim.res_ins_num ++;
              index = rand() % 8;
              if (index >= 4) {
                mut.read_seq[read_offset] = nt;
              } else {
                mut.read_seq[read_offset] = mut.ins_nt[index];
              }
              mut.maf_seq[maf_offset] = mut.read_seq[read_offset];
              mut.maf_ref_seq[maf_offset] = '-';
              read_offset ++;
            } else {
              mut.err_num ++;
              sim.res_del_num ++;
              mut.maf_seq[maf_offset] = '-';
              mut.maf_ref_seq[maf_offset] = nt;
              ref_offset ++;
            }
            maf_offset ++;
          }
        }
        mut.read_seq[read_offset] = '\0';
        mut.maf_seq[maf_offset] = '\0';
        mut.maf_ref_seq[maf_offset] = '\0';

        len = strlen(mut.read_seq);
        sim.res_len_total += len;
        freq_len[len] ++;

        if (len > sim.res_len_max) {
          sim.res_len_max = len;
        }
        if (len < sim.res_len_min) {
          sim.res_len_min = len;
        }

        value = 1.0 - ((double)mut.err_num / len);
        accuracy_total += value;
        accuracy = (int)(value * 100000 + 0.5);
        freq_accuracy[accuracy] ++;

        for (i=0; i<len; i++) {
          mut.new_qc[i] = '!';
        }
        mut.new_qc[len] = '\0';

        if (sim.pass_num == 1) {
          sprintf(id, "%s_%ld", sim.id_prefix, sim.res_num);
          fprintf(fp_fq, "@%s\n%s\n+%s\n%s\n", id, mut.read_seq, id, mut.new_qc);
        } else {
          sprintf(id, "%s/%ld/%ld", sim.id_prefix, sim.res_num, h);
          fprintf(fp_sam, "%s\t4\t*\t0\t255\t*\t*\t0\t0\t%s\t%s", id, mut.read_seq, mut.new_qc);
          fprintf(fp_sam, "\tcx:i:3\tip:B:C");
          for (i=0; i<len; i++) {
            fprintf(fp_sam, ",9");
          }
          fprintf(fp_sam, "\tnp:i:1\tpw:B:C");
          for (i=0; i<len; i++) {
            fprintf(fp_sam, ",9");
          }
          qeval = len - 1;
          fprintf(fp_sam, "\tqs:i:0\tqe:i:%ld\trq:f:%f\tsn:B:f,10.0,10.0,10.0,10.0\tzm:i:%ld\tRG:Z:ffffffff\n", qeval, sim.accuracy_mean, sim.res_num);
        }

        digit_num1[0] = 3;
        digit_num2[0] = 1 + count_digit(sim.res_num);
        digit_num[0] = (digit_num1[0] >= digit_num2[0]) ? digit_num1[0] : digit_num2[0];

        digit_num1[1] = count_digit((mut.seq_left - 1));
        digit_num2[1] = 1;
        digit_num[1] = (digit_num1[1] >= digit_num2[1]) ? digit_num1[1] : digit_num2[1];

        digit_num1[2] = count_digit((mut.seq_right - mut.seq_left + 1));
        digit_num2[2] = count_digit(len);
        digit_num[2] = (digit_num1[2] >= digit_num2[2]) ? digit_num1[2] : digit_num2[2];

        digit_num1[3] = count_digit(templ.len);
        digit_num2[3] = count_digit(len);
        digit_num[3] = (digit_num1[3] >= digit_num2[3]) ? digit_num1[3] : digit_num2[3];

        fprintf(fp_maf, "a\ns %s", templ.id);
        while (digit_num1[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num1[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld", mut.seq_left - 1);
        while (digit_num1[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld +", mut.seq_right - mut.seq_left + 1);
        while (digit_num1[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n", templ.len, mut.maf_ref_seq);
        fprintf(fp_maf, "s %s", id);
        while (digit_num2[0] ++ < digit_num[0]) {
          fprintf(fp_maf, " ");
        }
        while (digit_num2[1] ++ < digit_num[1]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %d", 0);
        while (digit_num2[2] ++ < digit_num[2]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %c", len, mut.seq_strand);
        while (digit_num2[3] ++ < digit_num[3]) {
          fprintf(fp_maf, " ");
        }
        fprintf(fp_maf, " %ld %s\n\n", len, mut.maf_seq);
      }
    }

    if (ret_pointer == NULL) {
      break;
    }

    ret = trim(line);

    if (line[0] == '>') {
      strncpy(templ.id, line+1, REF_ID_LEN_MAX);
      templ.id[REF_ID_LEN_MAX] = '\0';
      offset = 0;
      while (ret != EXISTS_LINE_FEED) {
        if (fgets(line, BUF_SIZE, fp) == NULL) {
          break;
        }
        ret = trim(line);
      }
    } else {
      copy_size = strlen(line);
      memcpy(templ.seq + offset, line, copy_size);
      offset += copy_size;
    }
  }
  fclose(fp);

  sim.res_pass_num = sim.res_num * sim.pass_num;
  sim.res_len_mean = (double)sim.res_len_total / sim.res_pass_num;
  sim.res_accuracy_mean = accuracy_total / sim.res_pass_num;

  if (sim.res_pass_num == 1) {
    sim.res_len_sd = 0.0;
    sim.res_accuracy_sd = 0.0;
  } else {
    variance = 0.0;
    for (i=0; i<=sim.len_max; i++) {
      if (freq_len[i] > 0) {
        variance += pow((sim.res_len_mean - i), 2) * freq_len[i];
      }
    }
    sim.res_len_sd = sqrt(variance / sim.res_pass_num);

    variance = 0.0;
    for (i=0; i<=100000; i++) {
      if (freq_accuracy[i] > 0) {
        variance += pow((sim.res_accuracy_mean - i * 0.00001), 2) * freq_accuracy[i];
      }
    }
    sim.res_accuracy_sd = sqrt(variance / sim.res_pass_num);
  }

  return SUCCEEDED;
}

/////////////////////////////////////////////////////////////
// Function: print_sim_param - Print simulation parameters //
/////////////////////////////////////////////////////////////

void print_sim_param() {
  fprintf(stderr, ":::: Simulation parameters :::\n\n");

  if (sim.strategy == STRATEGY_WGS) {
    fprintf(stderr, "strategy : wgs\n");
  } else if (sim.strategy == STRATEGY_TRANS) {
    fprintf(stderr, "strategy : trans\n");
  } else {
    fprintf(stderr, "strategy : templ\n");
  }
  if (sim.method == METHOD_QS) {
    fprintf(stderr, "method : qshmm\n");
    fprintf(stderr, "qshmm : %s\n", qshmm.file);
  } else if (sim.method == METHOD_ERR) {
    fprintf(stderr, "method : errhmm\n");
    fprintf(stderr, "errhmm : %s\n", errhmm.file);
  } else {
    fprintf(stderr, "method : sample\n");
  }

  if (sim.strategy == STRATEGY_WGS) {
    fprintf(stderr, "genome : %s\n", genome.file);
  } else if (sim.strategy == STRATEGY_TRANS) {
    fprintf(stderr, "transcript : %s\n", transcript.file);
  } else {
    fprintf(stderr, "template : %s\n", templ.file);
  }

  fprintf(stderr, "prefix : %s\n", sim.prefix);
  fprintf(stderr, "id-prefix : %s\n", sim.id_prefix);

  if (sim.strategy == STRATEGY_WGS) {
    fprintf(stderr, "depth : %lf\n", sim.depth);
  }

  if (sim.strategy == STRATEGY_TEMPL) {
  } else if ((sim.method == METHOD_QS) || (sim.method == METHOD_ERR)) {
    fprintf(stderr, "length-mean : %f\n", sim.len_mean);
    fprintf(stderr, "length-sd : %f\n", sim.len_sd);
    fprintf(stderr, "length-min : %ld\n", sim.len_min);
    fprintf(stderr, "length-max : %ld\n", sim.len_max);
  } else {
    fprintf(stderr, "length-mean : (sample FASTQ)\n");
    fprintf(stderr, "length-sd : (sample FASTQ)\n");
    fprintf(stderr, "length-min : %ld\n", sim.len_min);
    fprintf(stderr, "length-max : %ld\n", sim.len_max);
  }

  if (sim.method != METHOD_ERR) {
    fprintf(stderr, "difference-ratio : %ld:%ld:%ld\n", sim.sub_ratio, sim.ins_ratio, sim.del_ratio);
  }

  fprintf(stderr, "seed : %d\n", sim.seed);

  if ((sim.method == METHOD_QS) || (sim.method == METHOD_ERR)) {
    fprintf(stderr, "accuracy-mean : %f\n", sim.accuracy_mean);
  } else {
    fprintf(stderr, "sample : %s\n", sample.file);
    fprintf(stderr, "sample-profile-id : %s\n", sim.profile_id);
    fprintf(stderr, "accuracy-mean : (sample FASTQ)\n");
    fprintf(stderr, "accuracy-sd : (sample FASTQ)\n");
    fprintf(stderr, "accuracy-min : %f\n", sim.accuracy_min);
    fprintf(stderr, "accuracy-max : %f\n", sim.accuracy_max);
  }

  fprintf(stderr, "pass_num : %d\n", sim.pass_num);
  fprintf(stderr, "hp-del-bias : %f\n", sim.hp_del_bias);
  fprintf(stderr, "\n");
}

/////////////////////////////////////////////////////////////////
// Function: set_mut - Set mutation parameters and varianeces  //
/////////////////////////////////////////////////////////////////

int set_mut() {
  int i;

  for (i=0; i<=93; i++) {
    mut.err_thre[i] = int(qc[i].prob * 1000000 + 0.5);
    mut.sub_thre[i] = int((qc[i].prob * sim.sub_rate) * 1000000 + 0.5);
    mut.ins_thre[i] = int((qc[i].prob * (sim.sub_rate + sim.ins_rate)) * 1000000 + 0.5);
    mut.del_thre[i] = int((qc[i].prob * sim.del_rate) / (1 + qc[i].prob * sim.del_rate) * 1000000 + 0.5);
  }

  mut.sub_nt_a = (char*)"TGC";
  mut.sub_nt_t = (char*)"AGC";
  mut.sub_nt_g = (char*)"ATC";
  mut.sub_nt_c = (char*)"ATG";
  mut.sub_nt_n = (char*)"ATGC";
  mut.ins_nt = (char*)"ATGC";

  if ((mut.err = (char *)malloc(sim.len_max + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.qc = (char *)malloc(sim.len_max + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.new_qc = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.tmp_qc = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.read_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.maf_ref_seq = (char *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  if ((mut.hp = (short *)malloc(sim.len_max * 2 + 1)) == 0) {
    fprintf(stderr, "ERROR: Cannot allocate memory.\n");
    return FAILED;
  }

  return SUCCEEDED;
}


////////////////////////////////////////////////////////////////
// Function: print_simulation_stats - Print Simulation Stats. //
////////////////////////////////////////////////////////////////

void print_simulation_stats() {
  if (sim.strategy == STRATEGY_WGS) {
    sim.res_depth = (double)sim.res_len_total / genome.len / sim.pass_num;
    fprintf(stderr, ":::: Simulation stats (ref.%ld) ::::\n\n", genome.num);
    fprintf(stderr, "read num. : %ld\n", sim.res_num);
    fprintf(stderr, "depth : %lf\n", sim.res_depth);
  } else {
    fprintf(stderr, ":::: Simulation stats ::::\n\n");
    fprintf(stderr, "read num. : %ld\n", sim.res_num);
  }
  fprintf(stderr, "read length mean (SD) : %f (%f)\n",
    sim.res_len_mean, sim.res_len_sd);
  fprintf(stderr, "read length min : %ld\n", sim.res_len_min);
  fprintf(stderr, "read length max : %ld\n", sim.res_len_max);
  fprintf(stderr, "read accuracy mean (SD) : %f (%f)\n",
    sim.res_accuracy_mean, sim.res_accuracy_sd);
  sim.res_sub_rate = (double)sim.res_sub_num / sim.res_len_total;
  sim.res_ins_rate = (double)sim.res_ins_num / sim.res_len_total;
  sim.res_del_rate = (double)sim.res_del_num / sim.res_len_total;
  fprintf(stderr, "substitution rate. : %f\n", sim.res_sub_rate);
  fprintf(stderr, "insertion rate. : %f\n", sim.res_ins_rate);
  fprintf(stderr, "deletion rate. : %f\n", sim.res_del_rate);
  fprintf(stderr, "\n");
}

///////////////////////////////////////////////////
// Function: set_qshmm - Set quality code model  //
///////////////////////////////////////////////////

int set_qshmm() {
  FILE *fp;
  char line[BUF_SIZE];
  char *tp;
  int accuracy, state, num;
  int i, j, k;

  if ((fp = fopen(qshmm.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", qshmm.file);
    return FAILED;
  }

  for (i=0; i<=ACCURACY_MAX; i++) {
    qshmm.exist_hmm[i] = 0;
    for (j=0; j<=STATE_MAX; j++) {
      qshmm.ip[i][j] = 0.0;
      for (k=0; k<=93; k++) {
        qshmm.ep[i][j][k] = 0.0;
      }
      for (k=0; k<=STATE_MAX; k++) {
        qshmm.tp[i][j][k] = 0.0;
      }
    }
  }

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    trim(line);
    tp = strtok(line, " ");
    accuracy = atoi(tp);
    qshmm.exist_hmm[accuracy] = 1;

    tp = strtok(NULL, " ");
    if (strcmp(tp, "IP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      tp = strtok(NULL, " ");
      qshmm.ip[accuracy][state] = atof(tp);

    } else if (strcmp(tp, "EP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        qshmm.ep[accuracy][state][num] = atof(tp);
        num ++;
        tp = strtok(NULL, " ");
      }

    } else if (strcmp(tp, "TP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        num ++;
        qshmm.tp[accuracy][state][num] = atof(tp);
        tp = strtok(NULL, " ");
      }
    }
  }
  fclose(fp);

  return SUCCEEDED;
}

//////////////////////////////////////////
// Function: set_errhmm - Set err model //
//////////////////////////////////////////

int set_errhmm() {
  FILE *fp;
  char line[BUF_SIZE];
  char *tp;
  int accuracy, state, num;
  int i, j, k;

  if ((fp = fopen(errhmm.file, "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open file: %s\n", errhmm.file);
    return FAILED;
  }

  for (i=0; i<=ACCURACY_MAX; i++) {
    errhmm.exist_hmm[i] = 0;
    errhmm.state_max[i] = 0;
    for (j=0; j<=STATE_MAX; j++) {
      errhmm.ip[i][j] = 0.0;
      for (k=0; k<=3; k++) {
        errhmm.ep[i][j][k] = 0.0;
      }
      for (k=0; k<=STATE_MAX; k++) {
        errhmm.tp[i][j][k] = 0.0;
      }
    }
  }
  errhmm.acc_min = 100;
  errhmm.acc_max = 0;

  while (fgets(line, BUF_SIZE, fp) != NULL) {
    trim(line);
    tp = strtok(line, " ");
    accuracy = atoi(tp);
    errhmm.exist_hmm[accuracy] = 1;
    if (errhmm.acc_min > accuracy) {
      errhmm.acc_min = accuracy;
    }
    if (errhmm.acc_max < accuracy) {
      errhmm.acc_max = accuracy;
    }

    tp = strtok(NULL, " ");
    if (strcmp(tp, "IP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      tp = strtok(NULL, " ");
      errhmm.ip[accuracy][state] = atof(tp);
      errhmm.state_max[accuracy] = state;

    } else if (strcmp(tp, "EP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        errhmm.ep[accuracy][state][num] = atof(tp);
        num ++;
        tp = strtok(NULL, " ");
      }

    } else if (strcmp(tp, "TP") == 0) {
      tp = strtok(NULL, " ");
      state = atoi(tp);
      num = 0;
      tp = strtok(NULL, " ");
      while (tp != NULL) {
        num ++;
        errhmm.tp[accuracy][state][num] = atof(tp);
        tp = strtok(NULL, " ");
      }
    }
  }
  fclose(fp);

  return SUCCEEDED;
}

///////////////////////////////////////////////////////
// Function: get_time_cpu - Get CPU time             //
///////////////////////////////////////////////////////

long get_time_cpu() {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  return ru.ru_utime.tv_sec;
}

///////////////////////////////////////////////////////
// Function: get_time - Get time                     //
///////////////////////////////////////////////////////

long get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec;
}

///////////////////////////////////////
// Function: print_help - Print help //
///////////////////////////////////////

void print_help() {
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: pbsim [options] \n\n");
  fprintf(stderr, " [general options]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "  --prefix             prefix of output files (sd).\n");
  fprintf(stderr, "  --id-prefix          prefix of read ID (S).\n");
  fprintf(stderr, "  --seed               for a pseudorandom number generator (Unix time).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for whole genome sequencing]\n\n");
  fprintf(stderr, "  --strategy           wgs\n");
  fprintf(stderr, "  --genome             FASTA format file (text file only).\n");
  fprintf(stderr, "  --depth              depth of coverage (20.0).\n");
  fprintf(stderr, "  --length-min         minimum length (100).\n");
  fprintf(stderr, "  --length-max         maximum length (1000000).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for transcriptome sequencing]\n\n");
  fprintf(stderr, "  --strategy           trans\n");
  fprintf(stderr, "  --transcript         original format file.\n");
  fprintf(stderr, "  --length-min         minimum length (100).\n");
  fprintf(stderr, "  --length-max         maximum length (1000000).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for template sequencing]\n\n");
  fprintf(stderr, "  --strategy           templ\n");
  fprintf(stderr, "  --template           FASTA format file (text file only).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for quality score model]\n\n");
  fprintf(stderr, "  --method             qshmm\n");
  fprintf(stderr, "  --qshmm              quality score model.\n");
  fprintf(stderr, "  --length-mean        mean length (9000.0).\n");
  fprintf(stderr, "  --length-sd          standard deviation of length (7000.0).\n");
  fprintf(stderr, "  --accuracy-mean      mean accuracy (0.85).\n");
  fprintf(stderr, "  --pass-num           number of sequencing passes (1).\n");
  fprintf(stderr, "  --difference-ratio   difference (error) ratio (6:55:39).\n");
  fprintf(stderr, "                       (substitution:insertion:deletion)\n");
  fprintf(stderr, "                       Each value must be 0-1000, e.g. 1000:1:0 is OK.\n");
  fprintf(stderr, "                       Note that the above default value is for PacBio RS II;\n");
  fprintf(stderr, "                       22:45:33 for PacBio Sequel and 39:24:36 for ONT are\n");
  fprintf(stderr, "                       recommended.\n");
  fprintf(stderr, "  --hp-del-bias        bias intensity of deletion in homopolymer (1).\n");
  fprintf(stderr, "                       The option specifies the deletion rate at 10-mer, where\n");
  fprintf(stderr, "                       the deletion rate at 1-mer is 1. The bias intensity from\n");
  fprintf(stderr, "                       1-mer to 10-mer is proportional to the length of the\n");
  fprintf(stderr, "                       homopolymer.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for error model]\n\n");
  fprintf(stderr, "  --method             errhmm\n");
  fprintf(stderr, "  --errhmm             error model.\n");
  fprintf(stderr, "  --length-mean        mean length (9000.0).\n");
  fprintf(stderr, "  --length-sd          standard deviation of length (7000.0).\n");
  fprintf(stderr, "  --accuracy-mean      mean accuracy (0.85).\n");
  fprintf(stderr, "  --pass-num           number of sequencing passes (1).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " [options for sample-based method]\n\n");
  fprintf(stderr, " Note that the method can only be used for wag strategy.\n\n");
  fprintf(stderr, "  --sample             FASTQ format file to sample (text file only).\n");
  fprintf(stderr, "  --sample-profile-id  sample (filtered) profile ID.\n");
  fprintf(stderr, "                       When using --sample, profile is stored;\n");
  fprintf(stderr, "                       'sample_profile_<ID>.fastq', and\n");
  fprintf(stderr, "                       'sample_profile_<ID>.stats' are created.\n");
  fprintf(stderr, "                       When not using --sample, profile is re-used.\n");
  fprintf(stderr, "                       Note that when profile is used, --length-min,max,\n");
  fprintf(stderr, "                       --accuracy-min,max would be the same as the profile.\n");
  fprintf(stderr, "  --accuracy-min       minimum accuracy (0.75).\n");
  fprintf(stderr, "  --accuracy-max       maximum accuracy (1.00).\n");
  fprintf(stderr, "  --difference-ratio   difference (error) ratio (6:55:39).\n");
  fprintf(stderr, "                       (substitution:insertion:deletion)\n");
  fprintf(stderr, "                       Each value must be 0-1000, e.g. 1000:1:0 is OK.\n");
  fprintf(stderr, "                       Note that the above default value is for PacBio RS II;\n");
  fprintf(stderr, "                       22:45:33 for PacBio Sequel and 39:24:36 for ONT are\n");
  fprintf(stderr, "                       recommended.\n");
  fprintf(stderr, "  --hp-del-bias        bias intensity of deletion in homopolymer (1).\n");
  fprintf(stderr, "                       The option specifies the deletion rate at 10-mer, where\n");
  fprintf(stderr, "                       the deletion rate at 1-mer is 1. The bias intensity from\n");
  fprintf(stderr, "                       1-mer to 10-mer is proportional to the length of the\n");
  fprintf(stderr, "                       homopolymer.\n");
  fprintf(stderr, "\n");
}

/////////////////////////////////////////
// Function: count_digit - count digit //
/////////////////////////////////////////

int count_digit(long num) {
  int digit = 1;
  int quotient;

  quotient = int(num / 10);

  while (quotient != 0) {
    digit ++;
    quotient = int(quotient / 10);
  }

  return digit;
}  

///////////////////////////////////////////////////////
// Function: revcomp - convert to reverse complement //
///////////////////////////////////////////////////////

void revcomp(char* str) {
  int i, len;
  char c;

  len = strlen(str);

  for(i=0; i<len/2; i++) {
    c = str[i];
    str[i] = str[len-i-1];
    str[len-i-1] = c;
  }

  for(i=0; i<len; i++) {
    if (str[i] == 'A') {
      str[i] = 'T';
    } else if (str[i] == 'T') {
      str[i] = 'A';
    } else if (str[i] == 'G') {
      str[i] = 'C';
    } else if (str[i] == 'C') {
      str[i] = 'G';
    }
  }
}

/////////////////////////////////////////////////
// Function: revshort - reverse array of short //
/////////////////////////////////////////////////

void revshort(short* str, long len) {
  int i;
  short s;

  for(i=0; i<len/2; i++) {
    s = str[i];
    str[i] = str[len-i-1];
    str[len-i-1] = s;
  }
}
