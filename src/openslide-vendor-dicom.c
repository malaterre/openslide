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
 *
 */

/*
 * DICOM support for VL Whole Slide Microscopy Image Storage (1.2.840.10008.5.1.4.1.1.77.1.6)
 *
 * quickhash comes from (0008,0018) SOP Instance UID
 *
 */

#include <config.h>

#include "openslide-private.h"
#include "openslide-decode-dicom.h"

#include <glib.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

//struct image {
//  int32_t fileno;
//  int32_t start_in_file;
//  int32_t length;
//  int32_t imageno;   // used only for cache lookup
//  int refcount;
//};
//
//struct tile {
//  struct image *image;
//
//  // location in the image
//  double src_x;
//  double src_y;
//};
struct dicom_wsmis_ops_data {
  //struct _openslide_dicomcache *tc;
  gchar **datafile_paths;
};

struct level {
  struct _openslide_level base;
  struct _openslide_dicom_level dicoml;
  struct _openslide_grid *grid;
};

static void destroy(openslide_t *osr) {
  struct dicom_wsmis_ops_data *data = osr->data;
#if 0
  _openslide_dicomcache_destroy(data->tc);
  g_slice_free(struct dicom_wsmis_ops_data, data);

  for (int32_t i = 0; i < osr->level_count; i++) {
    struct level *l = (struct level *) osr->levels[i];
    _openslide_grid_destroy(l->grid);
    g_slice_free(struct level, l);
  }
  g_free(osr->levels);
#endif
}

static bool read_tile(openslide_t *osr,
                      cairo_t *cr,
                      struct _openslide_level *level,
                      int64_t tile_col, int64_t tile_row,
                      void *arg,
                      GError **err) {
#if 0
  struct level *l = (struct level *) level;
  struct _openslide_dicom_level *dicoml = &l->dicoml;
  dicom *dicom = arg;

  // tile size
  int64_t tw = dicoml->tile_w;
  int64_t th = dicoml->tile_h;

  // cache
  struct _openslide_cache_entry *cache_entry;
  uint32_t *tiledata = _openslide_cache_get(osr->cache,
                                            level, tile_col, tile_row,
                                            &cache_entry);
  if (!tiledata) {
    tiledata = g_slice_alloc(tw * th * 4);
    if (!_openslide_dicom_read_tile(dicoml, dicom,
                                   tiledata, tile_col, tile_row,
                                   err)) {
      g_slice_free1(tw * th * 4, tiledata);
      return false;
    }

    // clip, if necessary
    if (!_openslide_dicom_clip_tile(dicoml, tiledata,
                                   tile_col, tile_row,
                                   err)) {
      g_slice_free1(tw * th * 4, tiledata);
      return false;
    }

    // put it in the cache
    _openslide_cache_put(osr->cache, level, tile_col, tile_row,
                         tiledata, tw * th * 4,
                         &cache_entry);
  }

  // draw it
  cairo_surface_t *surface = cairo_image_surface_create_for_data((unsigned char *) tiledata,
                                                                 CAIRO_FORMAT_ARGB32,
                                                                 tw, th,
                                                                 tw * 4);
  cairo_set_source_surface(cr, surface, 0, 0);
  cairo_surface_destroy(surface);
  cairo_paint(cr);

  // done with the cache entry, release it
  _openslide_cache_entry_unref(cache_entry);
#endif

  return true;
}

static bool paint_region(openslide_t *osr, cairo_t *cr,
                         int64_t x, int64_t y,
                         struct _openslide_level *level,
                         int32_t w, int32_t h,
                         GError **err) {
#if 0
  struct dicom_wsmis_ops_data *data = osr->data;
  struct level *l = (struct level *) level;

  dicom *dicom = _openslide_dicomcache_get(data->tc, err);
  if (dicom == NULL) {
    return false;
  }

  bool success = _openslide_grid_paint_region(l->grid, cr, dicom,
                                              x / l->base.downsample,
                                              y / l->base.downsample,
                                              level, w, h,
                                              err);
  _openslide_dicomcache_put(data->tc, dicom);

  return success;
#endif
}

static const struct _openslide_ops dicom_wsmis_ops = {
  .paint_region = paint_region,
  .destroy = destroy,
};

bool _openslide_dicom_is_wsmis(struct _openslide_dicom_wsmis *tl,
                                  int64_t dir) {
  return true;
}

static bool dicom_wsmis_detect(const char *filename,
                                struct _openslide_tifflike *tl,
                                GError **err) {
  // reject TIFFs
  if (tl) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "Is a TIFF file");
    return false;
  }

  // verify filename
  if (!g_str_has_suffix(filename, "DICOMDIR")) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "File does not have %s extension", "DICOMDIR");
    return false;
  }

  // verify existence
  if (!g_file_test(filename, G_FILE_TEST_EXISTS)) {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                "File does not exist");
    return false;
  }

  // ensure DICOM is wsmis
//  if (!_openslide_dicom_is_wsmis(tl, 0)) {
//    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
//                "DICOM is not wsmis");
//    return false;
//  }

  return true;
}

static bool dicom_wsmis_open(openslide_t *osr, const char *filename,
                              struct _openslide_tifflike *tl G_GNUC_UNUSED,
                              struct _openslide_hash *quickhash1, GError **err) {
  char *dirname = NULL;
  // get directory from filename
  dirname = g_strndup(filename, strlen(filename) - strlen("DICOMDIR"));

  struct _openslide_dicom * instance = _openslide_dicom_create(filename, err);

  struct dicom_wsmis_ops_data *data = g_slice_new0(struct dicom_wsmis_ops_data);
  osr->data = data;

  // set ops
  osr->ops = &dicom_wsmis_ops;

  char **datafile_paths = NULL;
  if(!_openslide_dicom_readindex(instance, dirname, &datafile_paths))
    {
    g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
      "Could not read DICOMDIR");
    return false;
    }
  _openslide_dicom_destroy(instance);

  // accumulate tiled levels
  char ** fullpath = datafile_paths;
  GPtrArray *level_array = g_ptr_array_new();
  while( *fullpath )
    {
    printf( "ICI: %s\n", *fullpath );

    // create level
    struct level *l = g_slice_new0(struct level);
    struct _openslide_dicom_level *dicoml = &l->dicoml;
    //if (!_openslide_tiff_level_init(tiff,
    //                                TIFFCurrentDirectory(tiff),
    //                                (struct _openslide_level *) l,
    //                                tiffl,
    //                                err)) {
    //  g_slice_free(struct level, l);
    //  goto FAIL;
    //}
    l->grid = _openslide_grid_create_simple(osr,
                                            dicoml->tiles_across,
                                            dicoml->tiles_down,
                                            dicoml->tile_w,
                                            dicoml->tile_h,
                                            read_tile);

    // add to array
    g_ptr_array_add(level_array, l);

    ++fullpath;
    }
  struct level **levels =
    (struct level **) g_ptr_array_free(level_array, false);

  osr->levels = (struct _openslide_level **) levels;

  return true;
}

const struct _openslide_format _openslide_format_dicom_wsmis = {
  .name = "dicom-wsmis",
  .vendor = "dicom-wsmis",
  .detect = dicom_wsmis_detect,
  .open = dicom_wsmis_open,
};
