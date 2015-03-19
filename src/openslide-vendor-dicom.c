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

#include <glib.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

struct dicom_wsmis_ops_data {
  struct _openslide_dicomcache *tc;
};

struct level {
  struct _openslide_level base;
//  struct _openslide_dicom_level dicoml;
  struct _openslide_grid *grid;
};

static void destroy(openslide_t *osr) {
#if 0
  struct dicom_wsmis_ops_data *data = osr->data;
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

static bool dicom_wsmis_open(openslide_t *osr,
                              const char *filename,
                              struct _openslide_dicomlike *tl,
                              struct _openslide_hash *quickhash1,
                              GError **err) {
  char *dirname = NULL;
  // get directory from filename
  dirname = g_strndup(filename, strlen(filename) - strlen("DICOMDIR"));

  struct _openslide_dicom * instance = _openslide_dicom_create(filename, err);
  if(!_openslide_dicom_readindex(instance, err))
    {
    return false;
    }
  _openslide_dicom_destroy(instance);

#if 0
  GPtrArray *level_array = g_ptr_array_new();

  // open dicom
  struct _openslide_dicomcache *tc = _openslide_dicomcache_create(filename);
  dicom *dicom = _openslide_dicomcache_get(tc, err);
  if (!dicom) {
    goto FAIL;
  }

  // accumulate tiled levels
  do {
    // confirm that this directory is tiled
    if (!dicomIsTiled(dicom)) {
      continue;
    }

    // confirm it is either the first image, or reduced-resolution
    if (dicomCurrentDirectory(dicom) != 0) {
      uint32_t subfiletype;
      if (!dicomGetField(dicom, dicomTAG_SUBFILETYPE, &subfiletype)) {
        continue;
      }

      if (!(subfiletype & FILETYPE_REDUCEDIMAGE)) {
        continue;
      }
    }

    // verify that we can read this compression (hard fail if not)
    uint16_t compression;
    if (!dicomGetField(dicom, dicomTAG_COMPRESSION, &compression)) {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "Can't read compression scheme");
      goto FAIL;
    };
    if (!dicomIsCODECConfigured(compression)) {
      g_set_error(err, OPENSLIDE_ERROR, OPENSLIDE_ERROR_FAILED,
                  "Unsupported dicom compression: %u", compression);
      goto FAIL;
    }

    // create level
    struct level *l = g_slice_new0(struct level);
    struct _openslide_dicom_level *dicoml = &l->dicoml;
    if (!_openslide_dicom_level_init(dicom,
                                    dicomCurrentDirectory(dicom),
                                    (struct _openslide_level *) l,
                                    dicoml,
                                    err)) {
      g_slice_free(struct level, l);
      goto FAIL;
    }
    l->grid = _openslide_grid_create_simple(osr,
                                            dicoml->tiles_across,
                                            dicoml->tiles_down,
                                            dicoml->tile_w,
                                            dicoml->tile_h,
                                            read_tile);

    // add to array
    g_ptr_array_add(level_array, l);
  } while (dicomReadDirectory(dicom));

  // sort tiled levels
  g_ptr_array_sort(level_array, width_compare);

  // set hash and properties
  struct level *top_level = level_array->pdata[level_array->len - 1];
  if (!_openslide_dicomlike_init_properties_and_hash(osr, tl, quickhash1,
                                                    top_level->dicoml.dir,
                                                    0,
                                                    err)) {
    goto FAIL;
  }

  // unwrap level array
  int32_t level_count = level_array->len;
  struct level **levels =
    (struct level **) g_ptr_array_free(level_array, false);
  level_array = NULL;

  // allocate private data
  struct dicom_wsmis_ops_data *data =
    g_slice_new0(struct dicom_wsmis_ops_data);

  // store osr data
  g_assert(osr->data == NULL);
  g_assert(osr->levels == NULL);
  osr->levels = (struct _openslide_level **) levels;
  osr->level_count = level_count;
  osr->data = data;
  osr->ops = &dicom_wsmis_ops;

  // put dicom handle and store dicomcache reference
  _openslide_dicomcache_put(tc, dicom);
  data->tc = tc;

  return true;

 FAIL:
  // free the level array
  if (level_array) {
    for (uint32_t n = 0; n < level_array->len; n++) {
      struct level *l = level_array->pdata[n];
      _openslide_grid_destroy(l->grid);
      g_slice_free(struct level, l);
    }
    g_ptr_array_free(level_array, true);
  }
  // free dicom
  _openslide_dicomcache_put(tc, dicom);
  _openslide_dicomcache_destroy(tc);
#endif
  return false;
}

const struct _openslide_format _openslide_format_dicom_wsmis = {
  .name = "dicom-wsmis",
  .vendor = "dicom-wsmis",
  .detect = dicom_wsmis_detect,
  .open = dicom_wsmis_open,
};
