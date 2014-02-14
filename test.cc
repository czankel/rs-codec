//
// Copyright (C) Chris Zankel. All rights reserved.
// This code is subject to U.S. and other copyright laws and
// intellectual property protections.
//
// The contents of this file are confidential and proprietary to Chris Zankel.
//
//

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

#include "rs-codec.h"

const int c_k = 48;
const int c_n = 56;
const int c_cluster_count = 100;
const int c_test_packet_size = 1400;

static const int c_line_count = 4;	// line repeats after c_line_count lines
static const int c_line_size = 128;	// k <= c_line_size

#include "stimulus.h"
#include "reference.h"

#define OLD


/* 
 * Dump the content of the local 'lines'.
 * Each byte of a line represents a single packet (first byte of the packet)
 */

static void dump_lines(const char *name,
    int count,
    int length,
    unsigned char *buffer)
{
  unsigned char *ptr = buffer;

  printf("unsigned char %s[%d][%d] = {\n", name, c_line_count, c_line_size);
  for (int line = 0; line < c_line_count; line++)
  {
    int i;
    puts("{");
    for (i = 0; i < c_line_size - 1; i++)
    {
      printf("%3d,", *ptr++);
      if (i % 32 == 31)
        putchar('\n');
    }
    printf("% 3d\n},\n", *ptr++);
  }
  puts("};\n");

}

/*
 * Create 'c_line_count' x 'c_line_size' number of random data an fill
 * the local 'lines' 
 */

static void create_stimulus()
{
  int seed;

  int fd = open("/dev/random", 0);
  if (fd < 0)
  {
    perror("open random");
    return;
  }
  if (read(fd, &seed, 4) != 4)
  {
    fprintf(stderr, "cannot set seed\n");
  }

  close(fd);

  srand(seed);

  unsigned char *buffer = new unsigned char[c_line_count * c_line_size];
  unsigned char *ptr = buffer;

  for (int i = 0; i < c_line_count * c_line_size; i++)
    *ptr++ = (unsigned char)((unsigned long long) rand()  * 255 / RAND_MAX);

  dump_lines("stimulus", c_line_count, c_line_size, buffer);

  delete[] buffer;
}

/* -------------------------------------------------------------------------- */

class cluster
{
  unsigned char *m_buffers;

  public:
  const int m_num_syspkt;
  const int m_num_parpkt;
  const int m_pktsize;
  const int m_data_offset;
  const int m_parity_offset;
  const int m_packet_size;
  const int m_buffer_size;

  // Create a buffer for 'num_syspkt' + 'num_parpkt' number of packets of 'pktsize' size and offset (for alignment tests)
  cluster(int num_syspkt, int num_parpkt, int pktsize, int data_offset, int parity_offset)
    : m_num_syspkt(num_syspkt)
      , m_num_parpkt(num_parpkt)
      , m_pktsize(pktsize)
      , m_data_offset(data_offset)
      , m_parity_offset(parity_offset)
      , m_packet_size(pktsize)
      , m_buffer_size((pktsize + (data_offset > parity_offset ? data_offset : parity_offset) + 3) & ~3)
      {
        m_buffers = new unsigned char[(num_syspkt+num_parpkt) * m_buffer_size];
      }

  ~cluster()
  {
    delete[] m_buffers;
  }

  unsigned char* get_buffer(int index, bool parity)
  {
    return &m_buffers[index * m_buffer_size + (parity ? m_parity_offset : m_data_offset) ];
  }

  // Fill the data buffers from the stimulus lines; copy the 'n-th' byte of each line to the n-th buffer, so that the
  // buffer has at maximum c_line_count different values that are repeated across the full buffer
  void fill_system_packets()
  {
    unsigned char (*line)[128] = stimulus;

    memset(m_buffers, 0, m_num_syspkt * m_buffer_size);

    for (int p = 0; p < m_num_syspkt; p++)
    {
      unsigned char *ptr = &m_buffers[p * m_buffer_size + m_data_offset];
      for (int l = 0; l < m_packet_size; l++)
      {
        *ptr++ = line[l % c_line_count][p];
      }
    }
  }

  // Clear the parity packets
  void clear_parity_packets()
  {
    for (int p = m_num_syspkt; p < m_num_parpkt; p++)
    {
      memset(m_buffers + p * m_buffer_size, 0, m_buffer_size);
    }
  }

  // Dump the first c_line_count lines from the packets
  void dump(const char *name)
  {
    unsigned char *buffer = new unsigned char[c_line_count * c_line_size];
    unsigned char *ptr = buffer;
    memset(buffer, 0, c_line_count * c_line_size);

    for (int l = 0; l < c_line_count; l++)
    {
      int i;
      for (i = 0; i < m_num_syspkt; i++)
        *ptr++ = m_buffers[i * m_buffer_size + l + m_data_offset];
      for (/* i = m_num_syspkt */; i < m_num_parpkt; i++)
        *ptr++ = m_buffers[i * m_buffer_size + l + m_parity_offset];
      ptr += c_line_size - i;
    }

    dump_lines(name, c_line_count, m_num_syspkt + m_num_parpkt, buffer);
    delete[] buffer;
  }

  // verify that the parity is correct
  bool verify(unsigned char (*line)[c_line_size], char *errline)
  {
    for (int l = 0; l < c_line_count; l++)
    {
      int i;
      for (i = 0; i < m_num_syspkt; i++)
      {
        if (line[l][i] != m_buffers[i * m_buffer_size + l + m_data_offset])
        {
          sprintf(errline,
              "Reference vector mismatch @system packet %d index %d!\n", i, l);
          return false;
        }
      }
      for (/* i = m_num_syspkt */; i < m_num_parpkt; i++)
      {
        if (line[l][i] != m_buffers[i * m_buffer_size + l + m_parity_offset])
        {
          sprintf(errline,
              "Reference vector mismatch @parity packet %d index %d!\n", i, l);
          return false;
        }
      }
    }
    return true;
  }

  bool compare(cluster *other, char *errline)
  {
    if (other->m_packet_size != m_packet_size)
    {
      fprintf(stderr, "ERROR: Data length mismatch!\n");
      return false;
    }
    if (m_num_syspkt != other->m_num_syspkt || m_num_parpkt != other->m_num_parpkt)
    {
      fprintf(stderr, "ERROR: Number of system or parity packets mismatch!\n");
      return false;
    }

    int i;
    for (i = 0; i < m_num_syspkt; i++)
    {
      unsigned char *p1 = m_buffers + i * m_buffer_size + m_data_offset;
      unsigned char *p2 =
        other->m_buffers + i * other->m_buffer_size + other->m_data_offset;
      for (int j = 0; j < m_packet_size; j++, p1++, p2++)
      {
        if (*p1 != *p2)
        {
          sprintf(errline, "System data mismatch @%d: %d != %d\n", j, *p1, *p2);
          return false;
        }
      }
    }
    for ( /* i */; i < m_num_syspkt + m_num_parpkt; i++)
    {
      unsigned char *p1 = m_buffers + i * m_buffer_size + m_parity_offset;
      unsigned char *p2 =
        other->m_buffers + i * other->m_buffer_size + other->m_parity_offset;
      for (int j = 0; j < m_packet_size; j++, p1++, p2++)
      {
        if (*p1 != *p2)
        {
          sprintf(errline, "Parity data mismatch @%d: %d != %d\n", j, *p1, *p2);
          return false;
        }
      }
    }
    return true;
  }
};


/* -------------------------------------------------------------------------- */

class Fec
{
  public:

    Fec(int k, int n) : m_k(k), m_n(n)
    {
      m_codec = rs_codec_create(k, n);
    }

    cluster *create_cluster(int packet_size, int data_offset, int parity_offset)
    {
      return new cluster(m_k, m_n-m_k, packet_size, data_offset, parity_offset);
    }

    void free_cluster(cluster *cl)
    {
      delete cl;
    }

#ifdef OLD
    void encode(cluster *cl)
    {
      unsigned char *buflist[cl->m_num_syspkt + cl->m_num_parpkt];
      for (int i = 0 ; i < cl->m_num_syspkt; i++)
      {
        buflist[i] = cl->get_buffer(i, /* parity: */ false);
      }
      for (int k = m_k; k < m_n; k++)
      {
        rs_codec_encode(
            m_codec,
            buflist,
            cl->get_buffer(k,/* parity: */ true),
            k,
            cl->m_pktsize);
      }
    }
#else
    void encode(cluster *cl)
    {
      cl->clear_parity_packets();
      unsigned char *parlist[cl->m_num_parpkt];
      for (int i = 0; i < cl->m_num_parpkt; i++)
      {
        parlist[i] = cl->get_buffer(i + cl->m_num_syspkt, /* parity: */ true);
      }

      for (int i = 0; i < cl->m_num_syspkt; i++)
      {
        rs_codec_encode_one(
            m_codec,
            parlist,
            cl->get_buffer(i, /* parity:*/ false),
            i,
            0,
            cl->m_pktsize);
      }
    }
#endif

  private:
    struct rs_codec_parms *m_codec;
    int m_k;
    int m_n;
};

void dump_ref(int k, int n)
{
  char name[80];

  sprintf(name, "fecref_%d_%d", k, n);

  Fec f(k, n);
  cluster *c = f.create_cluster(c_line_count, 0, 0); 
  c->fill_system_packets();
  f.encode(c);
  c->dump(name);
  delete c;
}


void speed_test(int data_offset, int parity_offset)
{
  Fec f(c_k, c_n);
  cluster *c[c_cluster_count];
  struct timespec ts_start, ts_end;
  unsigned long t;

  for (int i = 0; i < c_cluster_count; i++)
  {
    c[i] = f.create_cluster(c_test_packet_size, data_offset, parity_offset);
    c[i]->fill_system_packets();
  }

  //clock_gettime(CLOCK_MONOTONIC, &ts_start);
  for (int i = 0; i < c_cluster_count; i++)
  {
    f.encode(c[i]);
  }
  //clock_gettime(CLOCK_MONOTONIC, &ts_end);

  t  = ((ts_end.tv_sec - ts_start.tv_sec) * 1000000000UL + ts_end.tv_nsec - ts_start.tv_nsec + 500000UL) / 1000000UL;

  printf("Encoding %d clusters with k = %d, n = %d took %ldms\n", c_cluster_count, c_k, c_n, t);

  for (int i = 0; i < c_cluster_count; i++)
  {
    delete c[i];
  }
}

int main(int argc, const char **argv)
{
  if (argc < 2)
  {
    printf("Use this test with the following commands\n"
        "verify   - Encode a sequence of packets and compare a subsample\n"
        "           of the packet with a reference\n"
        "speed    - Run speed test across different aligned buffers\n");
    return 0;
  }

  // use 'init' to create the stimulus; replace the lines in the testvalues.h file
  if (argc >=2 && strncmp(argv[1], "init", 4) == 0)
  {
    create_stimulus();
    return 0;
  }

  // use 'ref' to create fec reference encoded lines
  if (argc >=2 && strncmp(argv[1], "ref", 3) == 0)
  {
    dump_ref(32, 38);
    dump_ref(48, 55);
    dump_ref(64, 72);
  }

  if (argc == 2 && strncmp(argv[1], "verify", 6) == 0)
  {
    char errline[200];

    Fec f1(32, 38);
    Fec f2(48, 55);
    Fec f3(64, 72);

    for (int size = 4; size < 60; size++)
    {
      // create base cluster to compare against
      cluster *c1ref = f1.create_cluster(size, 0, 0);
      cluster *c2ref = f2.create_cluster(size, 0, 0);
      cluster *c3ref = f3.create_cluster(size, 0, 0);
      c1ref->fill_system_packets();
      c2ref->fill_system_packets();
      c3ref->fill_system_packets();
      f1.encode(c1ref);
      f2.encode(c2ref);
      f3.encode(c3ref);

      for (int s = 0; s < 4; s++)
      {
        for (int p = 0; p < 4; p++)
        {
          // create clusters for different offsets
          cluster *c1 = f1.create_cluster(size, s, p);
          cluster *c2 = f2.create_cluster(size, s, p);
          cluster *c3 = f3.create_cluster(size, s, p);
          c1->fill_system_packets();
          c2->fill_system_packets();
          c3->fill_system_packets();
          f1.encode(c1);
          f2.encode(c2);
          f3.encode(c3);

          // compare with reference line
          errline[0] = 0;
          printf("Testing (%d,%d) size %d data-offset %d parity-offset %d ... %s\n%s",
              32, 38, size, s, p,
              (c1->verify(fecref_32_38, errline) && c1->compare(c1ref, errline)) ?
              "PASSED" : "FAILED", errline);
          errline[0] = 0;
          printf("Testing (%d,%d) size %d data-offset %d parity-offset %d ... %s\n%s",
              48, 55, size, s, p,
              (c2->verify(fecref_48_55, errline) && c2->compare(c2ref, errline)) ?
              "PASSED" : "FAILED", errline);
          errline[0] = 0;
          printf("Testing (%d,%d) size %d data-offset %d parity-offset %d ... %s\n%s",
              64, 72, size, s, p,
              (c3->verify(fecref_64_72, errline) && c3->compare(c3ref, errline)) ?
              "PASSED" : "FAILED", errline);
        }
      }
      delete c1ref;
      delete c2ref;
      delete c3ref;
    }
  }

  // speed test
  if (argc == 2 && strncmp(argv[1], "speed", 5) == 0)
  {
    for (int s = 0; s < 4; s++)
    {
      for (int p = 0; p < 4; p++)
      {
        printf("Offset: syspacket %d parpacket %d\n", s, p);
        speed_test(s, p);
      }
    }

    printf("Comparison:\n"
        "ARM i.MX27 @400MHz:   Encoding 100 clusters with k,n = 48,56: 1561ms\n"
        "ARM i.MX27 @400MHz:   Encoding 100 clusters with k,n = 48,56:  770ms (optimized)\n"
        "Intel Q6600 @2.40GHz: Encoding 100 clusters with k,n = 48,56:  70ms\n");
  }

  return 0;
}
