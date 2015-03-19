/*
 *  OpenSlide, a library for reading whole slide image files
 *
 *  Copyright (c) 2015 Mathieu Malaterre
 *  All rights reserved.
 *
 *  OpenSlide is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, version 2.1.
 *
 *  OpenSlide is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with OpenSlide. If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#include <config.h>

#include "openslide-private.h"
#include "openslide-decode-dicom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#include <byteswap.h>
#endif

/*
 * Implementation details: this is a dumb DICOM parser.
 * It has not built-in DICOM dictionary. This is a very limited implementation
 * of a Part-10 conforming implementation, which can only deals with WSMIS
 * instance and DICOMDIR index (it will not handle implicit and Big Endian
 * Explicit TS by design, which makes it a non-standard DICOM parser, but allow
 * a very simple implementation for OpenSlide).
 * Since it does not allow Implicit TS, it will be incapable of dealing with
 * Undefined Length VR:UN attributes which can only ever appears after an
 * Implicit -> Explicit conversion. This is technically impossible for WSMIS
 * instances (and thus not handled here).
 * 
 * Optimisations:
 * It will always parse everything, while some Defined Length Item / Sequence
 * could have been skipped. Since OpenSlide assumes direct file access, this
 * may be of little use.
 */
typedef struct
{
  FILE * stream;
  size_t max_len;
  size_t cur_pos;
} source;

void init_source( source *s, FILE * stream, uint32_t len )
{
  s->stream = stream;
  assert( len != (uint32_t)-1 );
  s->max_len = len;
  s->cur_pos = 0;
}

static size_t source_size(source * s)
{
  return s->max_len;
}
static size_t source_tell(source * s)
{
  return s->cur_pos;
}
static inline size_t min( size_t a, size_t b )
{
  return a < b ? a : b;
}
static bool source_read( source * s, char * out, size_t len )
{
  assert( s->cur_pos <= s->max_len );
  size_t llen = min( len , s->max_len - s->cur_pos );
  const size_t ret = fread( out, 1, len, s->stream );
  const bool b = ret == len;
  s->cur_pos += llen;
  assert( s->cur_pos <= s->max_len );
  return b && len == llen;
}
/* will skip up to len bytes (never more than max size) */
static bool source_skip( source * s, uint32_t len )
{
  assert( s->cur_pos <= s->max_len );
  const size_t llen = min( len , s->max_len - s->cur_pos );
  bool b = true;
  if( llen )
    {
    int ret = fseeko(s->stream, llen, SEEK_CUR );
    b = ret == 0;
    assert( b );
    s->cur_pos += llen;
    }
  assert( s->cur_pos <= s->max_len );
  return b /*&& len == llen*/;
}

static bool read_preamble(FILE * stream)
{
  int ret = fseeko(stream, 128, SEEK_SET );
  assert( ret == 0 );
  char buf[4];
  size_t n = fread( buf, sizeof *buf, sizeof buf, stream );
  assert( n == 4 );
  assert( strncmp( buf, "DICM", 4 ) == 0 );
  return true;
}

typedef uint32_t tag_t;
typedef uint16_t vr_t;
typedef uint32_t vl_t;

typedef union { uint16_t tags[2]; tag_t tag; } utag_t;
typedef union { char str[2]; vr_t vr; } uvr_t;
typedef union { char bytes[4]; vl_t vl; } uvl_t;
typedef union { char bytes[2]; uint16_t vl16; } uvl16_t;

static inline uint16_t get_group( tag_t tag )
{
  return (uint16_t)(tag >> 16);
}
static inline uint16_t get_element( tag_t tag )
{
  return (uint16_t)(tag & (uint16_t)0xffff);
}

#define MAKE_TAG(group,element) (group << 16 | element)

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define SWAP_TAG(t) t.tag = MAKE_TAG(t.tags[0], t.tags[1])
#else
#define SWAP_TAG(t) t.tags[0] = bswap_16( t.tags[0] ); t.tags[1] = bswap_16( t.tags[1] )
#endif

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define MAKE_VR(left,right) (right << 8 | left)
#else
#define MAKE_VR(left,right) (left << 8 | right)
#endif

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define SWAP_VL(vl)
#define SWAP_VL16(vl)
#else
#define SWAP_VL(vl) vl = bswap_32(vl)
#define SWAP_VL16(vl) vl = bswap_16(vl)
#endif

enum {
  E_INVALID = 0, /* Item, Item Delimitation Item & Sequence Delimitation Item */
  E_AE = MAKE_VR('A','E'),
  E_AS = MAKE_VR('A','S'),
  E_AT = MAKE_VR('A','T'),
  E_CS = MAKE_VR('C','S'),
  E_DA = MAKE_VR('D','A'),
  E_DS = MAKE_VR('D','S'),
  E_DT = MAKE_VR('D','T'),
  E_FL = MAKE_VR('F','L'),
  E_FD = MAKE_VR('F','D'),
  E_IS = MAKE_VR('I','S'),
  E_LO = MAKE_VR('L','O'),
  E_LT = MAKE_VR('L','T'),
  E_OB = MAKE_VR('O','B'),
  E_OD = MAKE_VR('O','D'),
  E_OF = MAKE_VR('O','F'),
  E_OW = MAKE_VR('O','W'),
  E_PN = MAKE_VR('P','N'),
  E_SH = MAKE_VR('S','H'),
  E_SL = MAKE_VR('S','L'),
  E_SQ = MAKE_VR('S','Q'),
  E_SS = MAKE_VR('S','S'),
  E_ST = MAKE_VR('S','T'),
  E_TM = MAKE_VR('T','M'),
  E_UI = MAKE_VR('U','I'),
  E_UL = MAKE_VR('U','L'),
  E_UN = MAKE_VR('U','N'),
  E_US = MAKE_VR('U','S'),
  E_UT = MAKE_VR('U','T'),
};

static inline bool isvr_valid( const uvr_t uvr )
{
  //assert( vr >= MAKE_VR('A','A') && vr <= MAKE_VR('Z','Z') );
  if( uvr.str[0] < 'A' || uvr.str[0] > 'Z' ) return false;
  if( uvr.str[1] < 'A' || uvr.str[1] > 'Z' ) return false;
  return true;
}

static inline bool isvr32( const vr_t vr )
{
  switch( vr )
    {
  case E_AE:
  case E_AS:
  case E_AT:
  case E_CS:
  case E_DA:
  case E_DS:
  case E_DT:
  case E_FD:
  case E_FL:
  case E_IS:
  case E_LO:
  case E_LT:
  case E_PN:
  case E_SH:
  case E_SL:
  case E_SS:
  case E_ST:
  case E_TM:
  case E_UI:
  case E_UL:
  case E_US:
    /* 16bits: */
    return false;
  case E_OB:
  case E_OD:
  case E_OF:
  case E_OW:
  case E_SQ:
  case E_UN:
  case E_UT:
    /* 32bits: */
    return true;
  default:
    /* parser error, or newer DICOM standard */
    /* return 32bits by default (required) */
    return true;
    }
}

typedef struct
{
  tag_t tag;
  vr_t vr;
  vl_t vl;
} data_element;

static inline bool tag_equal_to( const data_element * de, tag_t tag )
{
  return de->tag == tag;
}

static inline bool tag_is_lower( const data_element * de, tag_t tag )
{
  return de->tag < tag;
}

static inline bool is_start( const data_element * de )
{
  static const tag_t start = MAKE_TAG( 0xfffe,0xe000 );
  const bool b = de->tag == start;
  // can be undefined or defined length
  return b;
}
static inline bool is_end_item( const data_element * de )
{
  static const tag_t end_item = MAKE_TAG( 0xfffe,0xe00d );
  const bool b = de->tag == end_item;
  if( b ) assert( de->vl == 0 );
  return b;
}
static inline bool is_end_sq( const data_element * de )
{
  static const tag_t end_sq = MAKE_TAG( 0xfffe,0xe0dd );
  const bool b = de->tag == end_sq;
  if( b ) assert( de->vl == 0 );
  return b;
}
static inline bool is_encapsulated_pixel_data( const data_element * de )
{
  static const tag_t pixel_data = MAKE_TAG( 0x7fe0,0x0010 );
  const bool is_pixel_data = tag_equal_to(de, pixel_data);
  if( is_pixel_data )
    {
    assert( de->vl == (uint32_t)-1 );
    assert( de->vr == E_OB || de->vr == E_OW );
    return true;
    }
  return false;
}
static inline bool is_undef_len( const data_element * de )
{
  const bool b = de->vl == (uint32_t)-1;
  if( b ) assert( de->vr == E_SQ || is_encapsulated_pixel_data( de ) || is_start(de) );
  return b;
}
static inline uint32_t compute_len( const data_element * de )
{
  assert( !is_undef_len( de ) );
  const bool is32 = isvr32( de->vr );
  if( is32 )
    {
    return 4 /* tag */ + 4 /* VR */ + 4 /* VL */ + de->vl /* VL */;
    }
  return 4 /* tag */ + 4 /* VR/VL */ + de->vl /* VL */;
}
static inline uint32_t compute_undef_len( const data_element * de, uint32_t len )
{
  assert( is_undef_len( de ) );
  assert( len != (uint32_t)-1 );
  return 4 /* tag */ + 4 /* VR */ + 4 /* VL */ + len;
}
typedef struct {
  tag_t *tags;
  int ntags;
  int size;
} tag_path;

tag_path * create_tag_path()
{
  tag_path *ptr = (tag_path*)malloc( sizeof(tag_path) );
#if 0
  memset(ptr, 0, sizeof(tag_path));
#else
  ptr->tags = malloc( 16 * sizeof(tag_t) );
  ptr->ntags = 0;
  ptr->size = 16;
#endif
  return ptr;
}

void destroy_path(tag_path *tp)
{
  free( tp->tags );
  free( tp );
}

tag_path * clear_path(tag_path *tp)
{
  tp->ntags = 0;
  return tp;
}

tag_path * push_tag( tag_path * tp, tag_t t )
{
  if( tp->ntags < tp->size )
    {
    tag_t * ptr = tp->tags + tp->ntags;
    *ptr = t;
    tp->ntags++;
    }
  else
    {
    assert(0);
    }
  return tp;
}

tag_t pop_tag( tag_path * tp )
{
  tp->ntags--;
  const tag_t * ptr = tp->tags + tp->ntags;
  return *ptr;
}

bool tag_path_equal_to( tag_path * tp1, tag_path * tp2 )
{
  if( tp1->ntags == tp2->ntags )
    {
    int t;
    for( t = 0; t < tp1->ntags; ++t )
      {
      const tag_t * ptr1 = tp1->tags + t;
      const tag_t * ptr2 = tp2->tags + t;
      if( *ptr1 != *ptr2 ) return false;
      }
    // ok
    return true;
    }
  return false;
}

bool tag_path_match( tag_path * tp1, tag_path * tp2 )
{
  if( tp1->ntags >= tp2->ntags )
    {
    int t;
    for( t = 0; t < tp2->ntags; ++t )
      {
      const tag_t * ptr1 = tp1->tags + t;
      const tag_t * ptr2 = tp2->tags + t;
      if( *ptr1 != *ptr2 ) return false;
      }
    // ok
    return true;
    }
  return false;
}

int tag_path_length( tag_path * tp )
{
  return tp->ntags;
}

void print_path( tag_path * tp )
{
  fprintf(stdout, "Path: " );
  int t;
  for( t = 0; t < tp->ntags; ++t )
    {
    if( t != 0 ) fprintf(stdout, ">" );
    const tag_t * ptr = tp->tags + t;
    fprintf(stdout, "%04x,%04x", get_group(*ptr), get_element(*ptr) );
    }
  fprintf(stdout, "\n" );
}

typedef struct {
  tag_t *tags;
  int size; /* max number of tags */
  int *ntags;
  int nsets;
} tag_path_set;

tag_path_set * create_tag_path_set()
{
  tag_path_set *ptr = (tag_path_set*)malloc( sizeof(tag_path_set) );
  ptr->tags = malloc( 512 * sizeof(tag_t) );
  ptr->size = 512;
  ptr->ntags = malloc( 16 * sizeof(int) );
  ptr->nsets = 0;
  return ptr;
}

void destroy_path_set( tag_path_set * tps )
{
  free( tps->ntags );
  free( tps->tags );
  free( tps );
}

bool find_tag_path( tag_path_set *tps, tag_path *tp )
{
  int s;
  int total = 0;
  tag_path ref;
  for( s = 0; s < tps->nsets; ++s )
    {
    tag_t * ptr = tps->tags + total;
    ref.tags = ptr;
    ref.ntags = ref.size = tps->ntags[s];
    
    const bool b = tag_path_equal_to( &ref, tp );
    if( b ) return true;
    total += tps->ntags[s];
    }
  return false;
}

bool match_tag_path( tag_path_set *tps, tag_path *tp )
{
#if 1
  return true;
#else
  int s;
  int total = 0;
  tag_path ref;
  for( s = 0; s < tps->nsets; ++s )
    {
    tag * ptr = tps->tags + total;
    ref.tags = ptr;
    ref.ntags = ref.size = tps->ntags[s];
    
    const bool b = tag_path_match( &ref, tp );
    if( b ) return true;
    total += tps->ntags[s];
    }
  return false;
#endif
}

void add_tag_path( tag_path_set *tps, tag_path *tp )
{
  int s;
  int total = 0;
  for( s = 0; s < tps->nsets; ++s )
    {
    total += tps->ntags[s];
    }
  assert( total <= tps->size );
  tag_t * ptr = tps->tags + total;
  for( s = 0; s < tp->ntags; ++s )
    {
    ptr[s] = tp->tags[s];
    }
  tps->ntags[ tps->nsets ] = tp->ntags;
  tps->nsets++;
  assert( find_tag_path( tps, tp ) );
}

typedef struct dataset dataset;
struct dataset {
  tag_path *cur_tp;
  tag_path_set *tps;
  void (*handle_attribute)(/*const*/ dataset * ds, const data_element * de, source * value);
  void *data;
};

tag_path * get_tag_path( dataset * ds )
{
  return ds->cur_tp;
}

static void print_indent(const dataset * ds)
{
  int i;
  int indent = tag_path_length( ds->cur_tp );
  for( i = 0; i < indent - 1; ++i )
    {
    printf( "  " );
    }
}

// explicit
static bool read_explicit( data_element * de, FILE * stream )
{
  utag_t t;
  uvr_t vr;
  uvl_t vl;

  // Tag
  size_t n = fread( t.tags, sizeof *t.tags, 2, stream );
  if( n != 2 ) return false;
  SWAP_TAG(t);
  assert( tag_is_lower( de, t.tag ) );

  // Value Representation
  n = fread( vr.str, sizeof *vr.str, 2, stream );
  /* a lot of VR are not valid (eg: non-ASCII), however the standard may add
   * them in a future edition, so only exclude the impossible ones:
   * - 0x0
   */
  if( n != 2 || !isvr_valid(vr) ) return false;
  const bool is32 = isvr32( vr.vr );

  // padding and/or 16bits VL
  uvl16_t vl16;
  n = fread( vl16.bytes, sizeof *vl16.bytes, 2, stream );
  if( n != 2 ) return false;

  // Value Length
  if( is32 )
    {
    /* padding must be set to zero */
    if( vl16.vl16 != 0 ) return false;

    n = fread( vl.bytes, 1, 4, stream );
    if( n != 4 ) return false;
    SWAP_VL(vl.vl);
    }
  else
    {
    SWAP_VL16(vl16.vl16);
    vl.vl = vl16.vl16;
    }
  de->tag = t.tag;
  de->vr  = vr.vr;
  de->vl  = vl.vl;
  return true;
}

static bool read_explicit_undef( data_element * de, FILE * stream )
{
  utag_t t;
  uvr_t vr;
  uvl_t vl;

  // Tag
  size_t n = fread( t.tags, sizeof *t.tags, 2, stream );
  if( n != 2 ) return false;
  SWAP_TAG(t);
  assert( tag_is_lower( de, t.tag ) );

  static const tag_t item_end = MAKE_TAG( 0xfffe,0xe00d );
  if( t.tag == item_end )
    {
    // Special case:
    n = fread( vl.bytes, sizeof *vl.bytes, 4, stream );
    if( n != 4 || vl.vl != 0 ) return false;
    vr.vr = E_INVALID;
    }
  else
    {
    // Value Representation
    assert( get_group(t.tag) != 0xfffe );
    n = fread( vr.str, sizeof *vr.str, 2, stream );
    if( n != 2 || !isvr_valid(vr) ) return false;
    const bool is32 = isvr32( vr.vr );

    // padding and/or 16bits VL
    uvl16_t vl16;
    n = fread( vl16.bytes, sizeof *vl16.bytes, 2, stream );
    if( n != 2 ) return false;

    // Value Length
    if( is32 )
      {
      /* padding must be set to zero */
      if( vl16.vl16 != 0 ) return false;

      n = fread( vl.bytes, sizeof *vl.bytes, 4, stream );
      if( n != 4 ) return false;
      SWAP_VL(vl.vl);
      }
    else
      {
      SWAP_VL16(vl16.vl16);
      vl.vl = vl16.vl16;
      }
    }
  de->tag = t.tag;
  de->vr  = vr.vr;
  de->vl  = vl.vl;
  return true;
}

// implicit
static bool read_implicit( data_element * de, FILE * stream )
{
  utag_t t;
  uvl_t vl;
  // Tag
  size_t n = fread( t.tags, sizeof *t.tags, 2, stream );
  if( n != 2 ) return false;
  SWAP_TAG(t);
  assert( tag_is_lower( de, t.tag ) );

  // Value Length (always 32bits)
  n = fread( vl.bytes, 1, 4, stream );
  if( n != 4 ) return false;
  SWAP_VL(vl.vl);

  de->tag = t.tag;
  de->vr  = E_INVALID;
  de->vl  = vl.vl;
  return true;
}

static inline bool read_ul( FILE * stream, uint32_t * value )
{
  if( fread( (char*)value, sizeof *value, 1, stream ) != 1 )
    return false;
  SWAP_VL(*value);
  return true;
}

static bool read_meta( FILE * stream )
{
  data_element de = { 0 /* tag */ };
  if( !read_explicit( &de, stream ) ) return false;
  assert( de.tag == MAKE_TAG(0x2,0x0) );
  assert( de.vr == E_UL );
  assert( de.vl == 4 );
  // file meta group length
  uint32_t gl;
  if( !read_ul( stream, &gl ) ) return false;
  // for now skip the meta header:
  fseeko( stream, gl, SEEK_CUR );
  return true;
}

static void process_attribute( dataset *ds, data_element *de, FILE * stream )
{
  assert( !is_start(de) && !is_end_item(de) && !is_end_sq(de) );
  // TODO: do not set group length, they are deprecated
  if( is_undef_len(de) )
    {
    if( ds->handle_attribute ) (ds->handle_attribute)(ds, de, NULL);
    }
  else
    {
    source s;
    init_source( &s, stream, de->vl );
    if( ds->handle_attribute ) (ds->handle_attribute)(ds, de, &s);

    size_t siz = source_size(&s);
    // go to end of value:
    bool b = source_skip( &s, siz );
    assert( b );
    }
}

static uint32_t read_sq_undef( dataset * ds, FILE * stream );
static uint32_t read_encapsulated_pixel_data( dataset * ds, FILE * stream );

/* read a single undefined length Item */
static uint32_t read_item_undef( dataset * ds, FILE * stream )
{
  data_element de = { 0 /* tag */ };
  uint32_t itemlen = 0;
  do
    {
    // carefully read an Item End or an Explicit Data Element:
    const bool b = read_explicit_undef(&de, stream);
    assert( b );
    if( is_end_item(&de) )
      {
      // end of Item
      itemlen += 4 /* tags */ + 4 /* VL */;
      break;
      }

    push_tag( ds->cur_tp, de.tag );
    process_attribute( ds, &de, stream );
    if( is_undef_len(&de) )
      {
      // Either undefined SQ or encapsulated Pixel Data:
      const bool is_encapsulated = is_encapsulated_pixel_data( &de );
      if( is_encapsulated )
        {
        const uint32_t epdlen = read_encapsulated_pixel_data( ds, stream );
        itemlen += compute_undef_len( &de, epdlen );
        }
      else
        {
        assert( de.vr == E_SQ ); // IVRLE
        const uint32_t seqlen = read_sq_undef( ds, stream );
        itemlen += compute_undef_len( &de, seqlen );
        }
      }
    else
      {
      itemlen += compute_len( &de );
      }
    pop_tag( ds->cur_tp );
    } while( 1 );

  return itemlen;
}

/* read a single defined length Item of length: itemlen */
static void read_item_def( dataset * ds, FILE * stream, const uint32_t itemlen )
{
  uint32_t curlen = 0;
  data_element de = { 0 /* tag */ };
  while( curlen != itemlen )
    {
    assert( curlen <= itemlen );
    // read attribute
    const bool b = read_explicit(&de, stream);
    assert( b );
    push_tag( ds->cur_tp, de.tag );
    process_attribute( ds, &de, stream );
    if( is_undef_len(&de) )
      {
      // If undefined SQ or encapsulated Pixel Data:
      const bool is_encapsulated = is_encapsulated_pixel_data( &de );
      if( is_encapsulated )
        {
        const uint32_t epdlen = read_encapsulated_pixel_data( ds, stream );
        curlen += compute_undef_len( &de, epdlen );
        }
      else
        {
        const uint32_t seqlen = read_sq_undef( ds, stream );
        curlen += compute_undef_len( &de, seqlen );
        }
      }
    else
      {
      curlen += compute_len( &de );
      }
    pop_tag( ds->cur_tp );
    }
}

/* read a single undefined length SQ */
/* return the actual length of the sequence */
static uint32_t read_sq_undef( dataset * ds, FILE * stream )
{
  data_element de;
  uint32_t seqlen = 0;
  do
    {
    // restart tag for each new Item:
    de.tag = 0;
    // Item start
    const bool b = read_implicit( &de, stream );
    assert( b );
    if( is_end_sq(&de) )
      {
      /* end of SQ */
      assert( de.vl == 0 );
      seqlen += 4 /* tag */ + 4 /* vl */;
      break;
      }
    assert( is_start(&de) );

    if( is_undef_len(&de) )
      {
      const uint32_t itemlen = read_item_undef( ds, stream );
      seqlen += 4 /* tag */ + 4 /* vl */ + itemlen;
      }
    else
      {
      seqlen += 4 /* tag */ + 4 /* vl */ + de.vl;
      tag_path * tp = get_tag_path( ds );
      const bool b = match_tag_path( ds->tps, tp );
      if( !b )
        {
        // skip over an entire Item
        fseeko(stream, de.vl, SEEK_CUR );
        }
      else
        {
        read_item_def( ds, stream, de.vl );
        }
      }
    } while( 1 );

  // post-condition
  assert( is_end_sq(&de) );

  return seqlen;
}

/* Nested Icon SQ  (with encapsulated Pixel Data) */
static uint32_t read_encapsulated_pixel_data( dataset * ds, FILE * stream )
{
  data_element de;
  uint32_t epdlen = 0;
  do
    {
    /* restart for each Item */
    de.tag = 0;
    /* Item start */
    const bool b= read_implicit(&de, stream);
    assert( b );
    /* end of SQ */
    epdlen += 4 /* tag */ + 4 /* vl */;
    if( is_end_sq(&de) ) break;
    assert( is_start(&de) );

    fseeko(stream, de.vl, SEEK_CUR );
    epdlen += de.vl;
    } while( 1 );

  // post-condition
  assert( is_end_sq(&de) );

  return epdlen;
}

/* read a single defined length SQ */
static void read_sq_def( dataset * ds, FILE * stream, const uint32_t seqlen )
{
  uint32_t curlen = 0;
  data_element de;
  while( curlen != seqlen )
    {
    // illegal:
    assert( curlen <= seqlen );
    // restart for each Item:
    de.tag = 0;
    // start
    const bool b = read_implicit( &de, stream );
    assert( b && is_start(&de) );

    if( is_undef_len(&de) )
      {
      const size_t itemlen = read_item_undef( ds, stream );
      curlen += 4 /* tag */ + 4 /* VL */ + itemlen;
      }
    else
      {
      curlen += 4 /* tag */ + 4 /* VL */ + de.vl;
      tag_path * tp = get_tag_path( ds );
      const bool b = match_tag_path( ds->tps, tp );
      if( !b )
        {
        // skip over an entire Item
        fseeko(stream, de.vl, SEEK_CUR );
        }
      else
        {
        read_item_def( ds, stream, de.vl );
        }
      }
    }
}

static void handle_attribute( /*const*/ dataset * ds, const data_element * de, source * s )
{
  //char **filenames = (char**)ds->data;
  GSList * list = (GList*)ds->data;
  (void)ds;
  (void)s;
  uvr_t u;
  u.vr = de->vr;
  assert( de->vr != E_INVALID );
  //print_indent(ds);
  //printf( "%04x,%04x %c%c %d\n", get_group(de->tag), get_element(de->tag), u.str[0], u.str[1], de->vl );
  const bool b = find_tag_path( ds->tps, ds->cur_tp );
  if( b )
    {
    print_path( ds->cur_tp );
    //assert( s );
    size_t len = source_size(s);
    char buf[128];
    assert( len < 127 );
    bool b = source_read( s, buf, len );
    assert(b);
    buf[len] = 0;
    char * p = buf;
    while( *p )
      {
      if( *p == '\\' ) *p = '/';
      ++p;
      }
    printf( "%s\n", buf );
    list = g_slist_append (list, strdup(buf));
    printf( "%d\n", g_slist_length (list) );
    }
  ds->data = list;
}

static bool read_dataset( dataset * ds, FILE * stream )
{
  static const tag_t pixel_data = MAKE_TAG( 0x7fe0,0x0010 );
  data_element de = { 0 /* tag */ };
  while( read_explicit( &de, stream ) && tag_is_lower( &de, pixel_data ) )
    {
    assert( get_group( de.tag ) != 0xfffe );
    assert( get_group( de.tag ) <= 0x7fe0 );
    push_tag( ds->cur_tp, de.tag );
    if( is_undef_len(&de) )
      {
      if( de.vr != E_SQ )
        {
        assert( de.vr == E_UN ); // IVRLE !
        assert(0);
        }
      process_attribute( ds, &de, stream );
      read_sq_undef(ds, stream);
      }
    else
      {
      // FIXME: technically a Sequence could be stored with VR:UN, this does
      // not make sense in the WSI world, thus it is not handled.
      if( de.vr == E_SQ )
        {
        tag_path * tp = get_tag_path( ds );
        const bool b = match_tag_path( ds->tps, tp );
        if( !b )
          {
          // skip over an entire SQ
          assert( de.vr == E_SQ );
          fseeko(stream, de.vl, SEEK_CUR );
          }
        else
          {
          read_sq_def( ds, stream, de.vl );
          }
        }
      else
        {
        process_attribute( ds, &de, stream );
        }
      }
    pop_tag( ds->cur_tp );
    }
  if( !feof( stream ) )
    {
    //assert( is_encapsulated_pixel_data(&de) );
    assert( de.vr == E_OB || de.vr == E_OW );
    assert( is_undef_len(&de) || !is_undef_len(&de) );
    return true;
    }
  return false;
}

static bool read_dataset2( dataset * ds, FILE * stream )
{
  data_element de = { 0 /* tag */ };
  while( read_explicit( &de, stream ) )
    {
    assert( get_group( de.tag ) != 0xfffe );
    assert( get_group( de.tag ) <= 0x7fe0 );
    push_tag( ds->cur_tp, de.tag );
    if( is_undef_len(&de) )
      {
      if( de.vr != E_SQ )
        {
        assert( de.vr == E_UN ); // IVRLE !
        assert(0);
        }
      process_attribute( ds, &de, stream );
      read_sq_undef(ds, stream);
      }
    else
      {
      // FIXME: technically a Sequence could be stored with VR:UN, this does
      // not make sense in the WSI world, thus it is not handled.
      if( de.vr == E_SQ )
        {
        tag_path * tp = get_tag_path( ds );
        const bool b = match_tag_path( ds->tps, tp );
        if( !b )
          {
          // skip over an entire SQ
          assert( de.vr == E_SQ );
          fseeko(stream, de.vl, SEEK_CUR );
          }
        else
          {
          read_sq_def( ds, stream, de.vl );
          }
        }
      else
        {
        process_attribute( ds, &de, stream );
        }
      }
    pop_tag( ds->cur_tp );
    }
  assert( feof( stream ) );
  return true;
}

#if 0
typedef struct wsmis
{
  FILE * stream;
  dataset ds;
} wsmis;

wsmis * create_wsmis_from_file(const char * filename)
{
  wsmis *instance = (wsmis*)malloc( sizeof *instance );

  //instance->stream = fopen( filename, "rb" );
  //FILE *f = _openslide_fopen(filename, "rb", err);

  instance->ds.cur_tp = create_tag_path();
  instance->ds.tps = create_tag_path_set();

  // 5200,9230>0048,021a
    {
    tag_path *tp = create_tag_path();
    const tag_t t0 = MAKE_TAG(0x5200,0x9230);
    const tag_t t1 = MAKE_TAG(0x0048,0x021a);
    const tag_t t2 = MAKE_TAG(0x0048,0x021f);
    push_tag(
      push_tag(
        push_tag( tp, t0 ), t1 ), t2);
    add_tag_path( instance->ds.tps, tp );
    destroy_path( tp );
    }
  instance->ds.handle_attribute = handle_attribute;

  return instance;
}

int read_header( wsmis * instance )
{
  if(  !read_preamble( instance->stream )
    || !read_meta( instance->stream )
    || !read_dataset( &instance->ds, instance->stream ) )
    {
    return 0;
    }

  return 1;
}

void delete_wsmis( wsmis * instance )
{
  destroy_path( instance->ds.cur_tp );
  destroy_path_set( instance->ds.tps );

  fclose( instance->stream );
  free( instance );
}
#endif

//struct node {
//  char * filename;
//  char * next;
//};
struct _openslide_dicom {
  //char *filename;
  FILE * stream;
  dataset ds;
//  char **filenames;
//  struct node * root;
  GSList *filenames;
};

struct dicom_patient {
  struct dicom_patient * next;
  struct dicom_study * study;
};

struct dicom_study {
  struct dicom_study * next;
  struct dicom_series * series;
};

struct dicom_series {
  struct dicom_series * next;
  struct dicom_image * image;
};

struct dicom_image {
  struct dicom_image * next;
};

struct _openslide_dicom *_openslide_dicom_create(const char *filename,
                                                       GError **err) {
  struct _openslide_dicom *instance = NULL;
  FILE *stream = _openslide_fopen(filename, "rb", err);

  // allocate struct
  instance = g_slice_new0(struct _openslide_dicom);
  instance->stream = stream;

  instance->ds.cur_tp = create_tag_path();
  instance->ds.tps = create_tag_path_set();

  return instance;
}

bool _openslide_dicom_readindex(struct _openslide_dicom *instance, GError **err)
{
  // (0004,1500) CS [CDCAB791\CDCAB791\7A474CCD\CDCAB790 ]         # 36,1-8 Referenced File ID
  // 0004,1220>0004,1500
    {
    tag_path *tp = create_tag_path();
    const tag_t t0 = MAKE_TAG(0x0004,0x1220);
    const tag_t t1 = MAKE_TAG(0x0004,0x1500);
    push_tag( push_tag( tp, t0 ), t1 );
    add_tag_path( instance->ds.tps, tp );
    destroy_path( tp );
    }

  instance->filenames = NULL;
  instance->ds.data = instance->filenames;
  instance->ds.handle_attribute = handle_attribute;
  if(  !read_preamble( instance->stream )
    || !read_meta( instance->stream )
    || !read_dataset2( &instance->ds, instance->stream ) )
    {
    assert(0);
    return false;
    }
  guint len = g_slist_length (instance->ds.data);
    printf ("DEBUG: %d\n", len );
  GSList * fn;
  for (fn = instance->ds.data; fn != NULL;
    fn = g_slist_next (fn))
    {
    char *s = (char*)fn->data;
    //GString * ss = fn->data;
    printf ("DEBUG: %s\n", s );
    }
  g_slist_free (instance->ds.data);

  return true;
}

void _openslide_dicom_destroy(struct _openslide_dicom *instance) {
  if (instance == NULL) {
    return;
  }
  destroy_path( instance->ds.cur_tp );
  destroy_path_set( instance->ds.tps );

  fclose( instance->stream );

  g_slice_free(struct _openslide_dicom, instance);
}


/*
(0048,0105) SQ (Sequence with undefined length)                   # u/l,1 Optical Path
    (fffe,e000) na (Item with undefined length) 
      (0022,0016) SQ (Sequence with undefined length)               # u/l,1 Illumination
      (fffe,e0dd)                                                                       
      (0022,0019) SQ (Sequence with undefined length)               # u/l,1 Lenses Code 
        (fffe,e000) na (Item with defined length)                                       
          (0008,0100) SH [A-00118 ]                                 # 8,1 Code Value    
          (0008,0102) SH [SRT ]                                     # 4,1 Coding Scheme 
          (0008,0104) LO [Slide overview lens ]                     # 20,1 Code Meaning 
      (fffe,e0dd)  
*/
