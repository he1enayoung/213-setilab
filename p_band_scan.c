#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <pthread.h>  // pthread api

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

typedef struct {
    signal *sig;
    int filter_order;
    int num_bands;
    double bandwidth;
    double *band_power;
    int start_band;
    int end_band;
} thread_args_t;

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads\n");
}

double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}

void *worker(void *arg){

    thread_args_t *args = (thread_args_t *)arg;
    double coefficients[args->filter_order + 1];

    for (int b = args->start_band; b < args->end_band; b++){
            // // Make the filter
        generate_band_pass(args->sig->Fs,
                        b * args->bandwidth + 0.0001, // keep within limits
                        (b + 1) * args->bandwidth - 0.0001,
                        args->filter_order,
                        coefficients);
        hamming_window(args->filter_order,coefficients);

        // Convolve
        convolve_and_compute_power(args->sig->num_samples,
                                args->sig->data,
                                args->filter_order,
                                coefficients,
                                &(args->band_power[b]));
        
    }

    return NULL;
}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub, int num_threads) {

  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  //double filter_coeffs[filter_order + 1];
  double band_power[num_bands];
  
  pthread_t threads[num_threads];
  thread_args_t args[num_threads];

  int bands_per_thread = num_bands / num_threads;

  for (int thread = 0; thread < num_threads; thread++) {

    //set all the arguments
    args[thread].sig = sig;
    args[thread].filter_order = filter_order;
    args[thread].num_bands = num_bands;
    args[thread].bandwidth = bandwidth;
    args[thread].band_power = band_power;
    args[thread].start_band = thread * bands_per_thread;
    if (thread == num_threads - 1) {
    args[thread].end_band = num_bands;
    } else {
    args[thread].end_band = (thread + 1) * bands_per_thread;
    }

    pthread_create(&threads[thread], NULL, worker, &args[thread]);

    // // Make the filter
    // generate_band_pass(sig->Fs,
    //                    band * bandwidth + 0.0001, // keep within limits
    //                    (band + 1) * bandwidth - 0.0001,
    //                    filter_order,
    //                    filter_coeffs);
    // hamming_window(filter_order,filter_coeffs);

    // // Convolve
    // convolve_and_compute_power(sig->num_samples,
    //                            sig->data,
    //                            filter_order,
    //                            filter_coeffs,
    //                            &(band_power[band]));

  }

    for (int t = 0; t < num_threads; t++) {
            pthread_join(threads[t], NULL);
        }

    // Timing end
    //double end = get_seconds();

  unsigned long long tend = get_cycle_count();
  double end = get_seconds();
  printf("Parallel processing time: %lf seconds\n", end - start);

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
User time        %lf seconds\n\
System time      %lf seconds\n\
Page faults      %ld\n\
Page swaps       %ld\n\
Blocks of I/O    %ld\n\
Signals caught   %ld\n\
Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 7) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);

  int num_threads = atoi(argv[6]);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_threads)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}

