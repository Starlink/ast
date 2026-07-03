#if !defined( YAML_BACKEND_H_INCLUDED )
#define YAML_BACKEND_H_INCLUDED
/*
*+
*  Name:
*     yaml_backend.h

*  Type:
*     C include file.

*  Purpose:
*     Abstraction layer over libyaml or libfyaml for use by yamlchan.c.

*  Description:
*     This header provides a unified interface to either libyaml or libfyaml
*     for use by yamlchan.c. Exactly one of the preprocessor macros YAML or
*     FYAML must be defined when this header is included.
*
*     For libyaml (YAML defined): most interface functions are simple
*     preprocessor defines that map directly to libyaml functions.
*
*     For libfyaml (FYAML defined): interface functions are static inline
*     implementations using libfyaml's API, including callback adapters for
*     I/O.
*
*     When neither YAML nor FYAML is defined, this header defines nothing.

*  Types defined:
*     AstYamlParser    Opaque parser handle.
*     AstYamlEmitter   Opaque emitter handle.
*     AstYamlEvent     Flat event struct; fields anchor/tag/value/etc. at top
*                      level regardless of event type.
*     AstYamlStyle     Enum of node/scalar style hints.
*     AstYamlReadCb    I/O read callback type.
*     AstYamlWriteCb   I/O write callback type.

*  Functions defined (as macros or static inline):
*     astYamlParserInitialize   Initialise parser.
*     astYamlParserDelete       Release parser resources.
*     astYamlParserSetInput     Register I/O read callback.
*     astYamlParserParse        Parse one event into an AstYamlEvent.
*     astYamlParserGetError     Primary error message and 0-based line/col.
*     astYamlParserGetContext   Context description and 0-based line/col (or NULL).
*     astYamlEventDelete        Release event resources.
*     astYamlEmitterInitialize  Initialise emitter.
*     astYamlEmitterDelete      Release emitter resources.
*     astYamlEmitterSetOutput   Register I/O write callback.
*     astYamlEmitterSetIndent   Set indentation increment.
*     astYamlEmitterGetError    Error message string from emitter.
*     astYamlEmitStreamStart    Emit a stream-start event.
*     astYamlEmitStreamEnd      Emit a stream-end event.
*     astYamlEmitDocumentStart  Emit a document-start event.
*     astYamlEmitDocumentEnd    Emit a document-end event.
*     astYamlEmitScalar         Emit a scalar event.
*     astYamlEmitMappingStart   Emit a mapping-start event.
*     astYamlEmitMappingEnd     Emit a mapping-end event.
*     astYamlEmitSequenceStart  Emit a sequence-start event.
*     astYamlEmitSequenceEnd    Emit a sequence-end event.

*  Copyright:
*     Copyright (C) 2026 Science and Technology Facilities Council.
*     All Rights Reserved.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.

*  Authors:
*     EMB: E. Madison Bray (STScI)

*  History:
*     19-JUN-2026 (EMB):
*        Original version.
*-
*/


/* I/O callback types used by both backends.  Signatures match libyaml's
   yaml_read_handler_t and yaml_write_handler_t but use standard C types. */
typedef int (*AstYamlReadCb)(void *, unsigned char *, size_t, size_t *);
typedef int (*AstYamlWriteCb)(void *, unsigned char *, size_t);

/* Style hint for emit functions and the event style field.
   Values cover both node styles (block/flow for mappings and sequences) and
   scalar styles (plain, single-quoted, double-quoted).

   This tries to cover the possibilities supported by both backends. */
typedef enum {
   AST_YAML_STYLE_ANY = 0,       /* no preference; let the emitter choose */
   AST_YAML_STYLE_BLOCK,         /* block style (mapping or sequence) */
   AST_YAML_STYLE_FLOW,          /* flow style (mapping or sequence) */
   AST_YAML_STYLE_PLAIN,         /* plain (unquoted) scalar */
   AST_YAML_STYLE_SINGLE_QUOTED, /* single-quoted scalar */
   AST_YAML_STYLE_DOUBLE_QUOTED  /* double-quoted scalar */
} AstYamlStyle;

#if defined( YAML )

#include <yaml.h>
#include <stdlib.h>
#include <string.h>

/* Wrap yaml_parser_t so we can cache the combined error message string. */
typedef struct {
   yaml_parser_t yp;
   char *errmsg;
} AstYamlParser;

/* Emitter type maps directly to libyaml. */
typedef yaml_emitter_t AstYamlEmitter;

/* Flat event struct.  For libyaml, the underlying yaml_event_t is
   stored inside the struct and its fields are mirrored into the flat
   public fields by astYamlParserParse. */
typedef struct {
   int type;
   const char *anchor;
   const char *tag;
   const char *value;
   size_t value_length;
   int style;
   int plain_implicit;
   int quoted_implicit;
   yaml_event_t ye;
} AstYamlEvent;

/* Parser lifecycle.  Initialise also zeroes the errmsg cache. */
static inline int astYamlParserInitialize( AstYamlParser *p ) {
   memset( p, 0, sizeof(*p) );
   return yaml_parser_initialize( &p->yp );
}
static inline void astYamlParserDelete( AstYamlParser *p ) {
   yaml_parser_delete( &p->yp );
   free( p->errmsg );
   p->errmsg = NULL;
}
#define astYamlParserSetInput(p, cb, d) \
   yaml_parser_set_input( &(p)->yp, (yaml_read_handler_t *)(cb), d )

/* Error detail accessors.  Line and column are 0-based; line is -1 when
   the location is not available.
   astYamlParserGetError returns a heap-allocated string that combines
   problem and (when present) context; the string is owned by the parser
   and freed on the next call or when the parser is deleted.
   astYamlParserGetContext always returns NULL for libyaml because context
   is already folded into the combined message. */
static inline const char *astYamlParserGetError( AstYamlParser *p,
                                                  int *line, int *col ) {
   const char *problem = p->yp.problem ? p->yp.problem : "unknown error";
   const char *context = p->yp.context;
   char *msg;
   size_t len;
   free( p->errmsg );
   p->errmsg = NULL;
   *line = (int) p->yp.problem_mark.line;
   *col  = (int) p->yp.problem_mark.column;
   if( context ) {
      len = strlen( problem ) + strlen( context ) + 4; /* " (" + ")" + NUL */
      msg = malloc( len );
      if( msg ) {
         snprintf( msg, len, "%s (%s)", problem, context );
         p->errmsg = msg;
         return msg;
      }
   }
   return problem;
}
static inline const char *astYamlParserGetContext( AstYamlParser *p,
                                                    int *line, int *col ) {
   (void) p;
   *line = -1;
   *col  = -1;
   return NULL;
}

/* Event delete: operate on the embedded yaml_event_t. */
#define astYamlEventDelete(ev) yaml_event_delete( &(ev)->ye )

/* Emitter lifecycle: map directly to libyaml, except astYamlEmitterInitialize
   also sets the line break style to Unix (\n) which is always the right
   choice. */
static inline int astYamlEmitterInitialize( AstYamlEmitter *e ) {
   if( !yaml_emitter_initialize( e ) )
      return 0;

   yaml_emitter_set_break( e, YAML_LN_BREAK );
   return 1;
}
#define astYamlEmitterDelete(e) yaml_emitter_delete( e )
#define astYamlEmitterSetOutput(e, cb, d) \
   yaml_emitter_set_output( e, (yaml_write_handler_t *)(cb), d )
#define astYamlEmitterSetIndent(e, n) yaml_emitter_set_indent( e, n )
static inline const char *astYamlEmitterGetError( AstYamlEmitter *e ) {
   return e->problem;
}

/* astYamlParserParse: parse one event and fill the flat AstYamlEvent. */
static inline int astYamlParserParse( AstYamlParser *parser,
                                      AstYamlEvent *ev ) {
   int ok;
   ok = yaml_parser_parse( &parser->yp, &ev->ye );
   ev->type = ok ? ev->ye.type : YAML_NO_EVENT;
   ev->anchor = NULL;
   ev->tag = NULL;
   ev->value = NULL;
   ev->value_length = 0;
   ev->style = 0;
   ev->plain_implicit = 0;
   ev->quoted_implicit = 0;

   if( ok ) {
      switch( ev->ye.type ) {
      case YAML_SCALAR_EVENT:
         ev->anchor = (const char *) ev->ye.data.scalar.anchor;
         ev->tag = (const char *) ev->ye.data.scalar.tag;
         ev->value = (const char *) ev->ye.data.scalar.value;
         ev->value_length = ev->ye.data.scalar.length;
         ev->style = ev->ye.data.scalar.style;
         ev->plain_implicit = ev->ye.data.scalar.plain_implicit;
         ev->quoted_implicit = ev->ye.data.scalar.quoted_implicit;
         break;
      case YAML_MAPPING_START_EVENT:
         ev->anchor = (const char *) ev->ye.data.mapping_start.anchor;
         ev->tag = (const char *) ev->ye.data.mapping_start.tag;
         ev->style = ev->ye.data.mapping_start.style;
         ev->plain_implicit = ev->ye.data.mapping_start.implicit;
         break;
      case YAML_SEQUENCE_START_EVENT:
         ev->anchor = (const char *) ev->ye.data.sequence_start.anchor;
         ev->tag = (const char *) ev->ye.data.sequence_start.tag;
         ev->style = ev->ye.data.sequence_start.style;
         ev->plain_implicit = ev->ye.data.sequence_start.implicit;
         break;
      case YAML_ALIAS_EVENT:
         ev->anchor = (const char *) ev->ye.data.alias.anchor;
         break;
      default:
         break;
      }
   }
   return ok;
}

/* Translate AstYamlStyle to libyaml scalar/node style enums. */
static inline yaml_scalar_style_t _astYamlToScalarStyle( AstYamlStyle s ) {
   switch( s ) {
   case AST_YAML_STYLE_PLAIN:
      return YAML_PLAIN_SCALAR_STYLE;
   case AST_YAML_STYLE_SINGLE_QUOTED:
      return YAML_SINGLE_QUOTED_SCALAR_STYLE;
   case AST_YAML_STYLE_DOUBLE_QUOTED:
      return YAML_DOUBLE_QUOTED_SCALAR_STYLE;
   default:
      return YAML_ANY_SCALAR_STYLE;
   }
}

static inline yaml_mapping_style_t _astYamlToMappingStyle( AstYamlStyle s ) {
   switch( s ) {
   case AST_YAML_STYLE_BLOCK:
      return YAML_BLOCK_MAPPING_STYLE;
   case AST_YAML_STYLE_FLOW:
      return YAML_FLOW_MAPPING_STYLE;
   default:
      return YAML_ANY_MAPPING_STYLE;
   }
}

static inline yaml_sequence_style_t _astYamlToSequenceStyle( AstYamlStyle s ) {
   switch( s ) {
   case AST_YAML_STYLE_BLOCK:
      return YAML_BLOCK_SEQUENCE_STYLE;
   case AST_YAML_STYLE_FLOW:
      return YAML_FLOW_SEQUENCE_STYLE;
   default:
      return YAML_ANY_SEQUENCE_STYLE;
   }
}

/* Emit functions: combine event initialisation and emission.
   Return 1 on success, 0 on failure. */

static inline int astYamlEmitStreamStart( AstYamlEmitter *emitter ) {
   yaml_event_t ev;
   yaml_stream_start_event_initialize( &ev, YAML_UTF8_ENCODING );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitStreamEnd( AstYamlEmitter *emitter ) {
   yaml_event_t ev;
   yaml_stream_end_event_initialize( &ev );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitDocumentStart( AstYamlEmitter *emitter,
                                            int with_version,
                                            int with_tags,
                                            const char *tag_handle,
                                            const char *tag_prefix ) {
   yaml_event_t ev;
   yaml_version_directive_t ver = { 1, 1 };
   yaml_tag_directive_t tag;

   if( with_tags ) {
      tag.handle = (yaml_char_t *) tag_handle;
      tag.prefix = (yaml_char_t *) tag_prefix;
   }

   yaml_document_start_event_initialize( &ev,
      with_version ? &ver : NULL,
      with_tags ? &tag : NULL,
      with_tags ? &tag + 1 : NULL,
      0 );

   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitDocumentEnd( AstYamlEmitter *emitter,
                                          int implicit ) {
   yaml_event_t ev;
   yaml_document_end_event_initialize( &ev, implicit );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitScalar( AstYamlEmitter *emitter,
                                     const char *anchor,
                                     const char *tag, const char *value,
                                     size_t len,
                                     int plain_implicit, int quoted_implicit,
                                     AstYamlStyle style ) {
   yaml_event_t ev;
   yaml_scalar_event_initialize( &ev,
      (yaml_char_t *) anchor, (yaml_char_t *) tag,
      (yaml_char_t *) value, len,
      plain_implicit, quoted_implicit, _astYamlToScalarStyle( style ) );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitMappingStart( AstYamlEmitter *emitter,
                                           const char *anchor,
                                           const char *tag, int implicit,
                                           AstYamlStyle style ) {
   yaml_event_t ev;
   yaml_mapping_start_event_initialize( &ev,
      (yaml_char_t *) anchor, (yaml_char_t *) tag,
      implicit, _astYamlToMappingStyle( style ) );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitMappingEnd( AstYamlEmitter *emitter ) {
   yaml_event_t ev;
   yaml_mapping_end_event_initialize( &ev );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitSequenceStart( AstYamlEmitter *emitter,
                                            const char *anchor,
                                            const char *tag, int implicit,
                                            AstYamlStyle style ) {
   yaml_event_t ev;
   yaml_sequence_start_event_initialize( &ev,
      (yaml_char_t *) anchor, (yaml_char_t *) tag,
      implicit, _astYamlToSequenceStyle( style ) );
   return yaml_emitter_emit( emitter, &ev );
}

static inline int astYamlEmitSequenceEnd( AstYamlEmitter *emitter ) {
   yaml_event_t ev;
   yaml_sequence_end_event_initialize( &ev );
   return yaml_emitter_emit( emitter, &ev );
}

#elif defined( FYAML )

#include <libfyaml.h>
#include <stdlib.h>
#include <string.h>

/* Event type constants: map YAML_* names to equivalent FYET_* values so that
   the switch/if chains in yamlchan.c need no changes. */
#define YAML_NO_EVENT             FYET_NONE
#define YAML_STREAM_START_EVENT   FYET_STREAM_START
#define YAML_STREAM_END_EVENT     FYET_STREAM_END
#define YAML_DOCUMENT_START_EVENT FYET_DOCUMENT_START
#define YAML_DOCUMENT_END_EVENT   FYET_DOCUMENT_END
#define YAML_ALIAS_EVENT          FYET_ALIAS
#define YAML_SCALAR_EVENT         FYET_SCALAR
#define YAML_SEQUENCE_START_EVENT FYET_SEQUENCE_START
#define YAML_SEQUENCE_END_EVENT   FYET_SEQUENCE_END
#define YAML_MAPPING_START_EVENT  FYET_MAPPING_START
#define YAML_MAPPING_END_EVENT    FYET_MAPPING_END

#define YAML_STR_TAG "tag:yaml.org,2002:str"

/* Parser handle: wraps struct fy_parser plus a buffered input copy
   that must outlive the parser (fy_parser_set_string does not copy). */
typedef struct {
   struct fy_parser *fyp;
   char *buf;
   AstYamlReadCb read_cb;
   void *read_data;
} AstYamlParser;

/* Emitter handle: wraps struct fy_emitter; the fy_emitter is created
   lazily at the first emit call so that indent can be set beforehand.
   line_buf accumulates partial output fragments between newlines so that
   AstYamlWriter always receives complete newline-terminated lines. */
typedef struct {
   struct fy_emitter *fye;
   AstYamlWriteCb write_cb;
   void *write_data;
   int indent;
   unsigned char *line_buf;
   size_t line_buf_len;
   size_t line_buf_cap;
} AstYamlEmitter;

/* Flat event struct.  For libfyaml the underlying struct fy_event * is
   stored so that it can be freed by astYamlEventDelete. */
typedef struct {
   int type;
   const char *anchor;
   const char *tag;
   const char *value;
   size_t value_length;
   int style;
   int plain_implicit;
   int quoted_implicit;
   struct fy_event *fye;
   struct fy_parser *fyp;
} AstYamlEvent;

/* Callback adapter: adapts AstYamlWriteCb to libfyaml's output fn type.
   libfyaml calls this with individual tokens rather than whole lines, but
   AstYamlWriter expects complete newline-terminated lines.  We therefore
   buffer fragments in emitter->line_buf and only forward to write_cb when a
   complete line (ending with '\n') has been assembled.  Returns the number of
   bytes consumed, or -1 on error. */
static int _astFyamlWriteAdapter( struct fy_emitter *emit,
                                  enum fy_emitter_write_type type,
                                  const char *str, int len, void *userdata ) {
   AstYamlEmitter *emitter = (AstYamlEmitter *) userdata;
   const unsigned char *p = (const unsigned char *) str;
   const unsigned char *end = p + (size_t) len;
   const unsigned char *nl;
   size_t n;
   unsigned char *new_buf;

   (void) emit; /* not used */
   (void) type; /* not used */

   while( p < end ) {
      nl = memchr( p, '\n', (size_t)( end - p ) );
      n  = nl ? (size_t)( nl - p + 1 ) : (size_t)( end - p );

      if( emitter->line_buf_len + n > emitter->line_buf_cap ) {
         size_t new_cap = emitter->line_buf_cap * 2 + n + 256;
         new_buf = realloc( emitter->line_buf, new_cap );
         if( !new_buf ) return -1;
         emitter->line_buf = new_buf;
         emitter->line_buf_cap = new_cap;
      }
      memcpy( emitter->line_buf + emitter->line_buf_len, p, n );
      emitter->line_buf_len += n;

      if( nl ) {
         /* Pass the line WITHOUT its trailing '\n'.  AstYamlWriter loops
            to pend and emits the content via astPutNextText, which adds
            its own '\n'.  If we included the '\n' here, AstYamlWriter
            would emit the content at the '\n' and then again at pend
            (as an empty string), producing a spurious blank line. */
         if( !emitter->write_cb( emitter->write_data,
                                 emitter->line_buf,
                                 emitter->line_buf_len - 1 ) )
            return -1;
         emitter->line_buf_len = 0;
      }

      p += n;
   }

   return len;
}

/* Error detail accessors using libfyaml's fy_diag error-collection API.
   Line and column come from the first collected error; line is -1 when
   not available.  Context is always NULL (libfyaml has no equivalent). */
static inline const char *astYamlParserGetError( AstYamlParser *p,
                                                  int *line, int *col ) {
   struct fy_diag *diag;
   struct fy_diag_error *err;
   void *iter = NULL;
   *line = -1; *col = -1;
   if( !p->fyp ) return "libfyaml parser not initialised";
   diag = fy_parser_get_diag( p->fyp );
   if( !diag ) return "unknown error";
   err = fy_diag_errors_iterate( diag, &iter );
   fy_diag_unref( diag );
   if( !err || !err->msg ) return "unknown error";
   /* fy_diag stores line/column 1-based; subtract 1 to match libyaml's
      0-based convention (callers add 1 back for display). */
   *line = err->line - 1;
   *col  = err->column - 1;
   return err->msg;
}
static inline const char *astYamlParserGetContext( AstYamlParser *p,
                                                    int *line, int *col ) {
   (void) p;
   *line = -1; *col = -1;
   return NULL;
}
static inline const char *astYamlEmitterGetError( AstYamlEmitter *e ) {
   struct fy_diag *diag;
   struct fy_diag_error *err;
   void *iter = NULL;
   if( !e->fye ) return "libfyaml emitter not initialised";
   diag = fy_emitter_get_diag( e->fye );
   if( !diag ) return "unknown error";
   err = fy_diag_errors_iterate( diag, &iter );
   fy_diag_unref( diag );
   if( !err || !err->msg ) return "unknown error";
   return err->msg;
}

/* Parser lifecycle. */

static inline int astYamlParserInitialize( AstYamlParser *parser ) {
   memset( parser, 0, sizeof( *parser ) );
   return 1;
}

static inline void astYamlParserDelete( AstYamlParser *parser ) {
   if( parser->fyp ) {
      fy_parser_destroy( parser->fyp );
      parser->fyp = NULL;
   }
   free( parser->buf );
   parser->buf = NULL;
}

/* astYamlParserSetInput: register a read callback.
   Because libfyaml's set-string API requires the buffer to remain valid
   for the lifetime of the parser (this allows it zero-copy reads, though
   we mostly use the copying *0 variants to get null-terminated strings),
   we buffer the entire input here. */
static inline void astYamlParserSetInput( AstYamlParser *parser,
                                          AstYamlReadCb cb,
                                          void *data ) {
   unsigned char tmp[ 4096 ];
   size_t nread;
   size_t buf_len = 0;
   size_t buf_cap = 0;
   char *new_buf;
   struct fy_parse_cfg cfg;

   parser->read_cb = cb;
   parser->read_data = data;

   while( 1 ) {
      nread = 0;
      if( !cb( data, tmp, sizeof( tmp ), &nread ) || nread == 0 ) {
         break;
      }

      if( buf_len + nread > buf_cap ) {
         buf_cap = buf_cap * 2 + nread + 1;
         new_buf = realloc( parser->buf, buf_cap );

         if( !new_buf ) {
            return;
         }

         parser->buf = new_buf;
      }

      memcpy( parser->buf + buf_len, tmp, nread );
      buf_len += nread;
   }

   memset( &cfg, 0, sizeof( cfg ) );
   parser->fyp = fy_parser_create( &cfg );

   if( parser->fyp ) {
      struct fy_diag *diag = fy_parser_get_diag( parser->fyp );
      if( diag ) {
         fy_diag_set_collect_errors( diag, true );
         fy_diag_unref( diag );
      }
      if( buf_len > 0 ) {
         fy_parser_set_string( parser->fyp, parser->buf, buf_len );
      }
   }
}

/* astYamlParserParse -- parse one event and fill the flat AstYamlEvent.
   NOTE: fy_token_get_text0 returns a pointer into the token; the string
   remains resident in memory until the event is freed by
   astYamlEventDelete. */
static inline int astYamlParserParse( AstYamlParser *parser,
                                      AstYamlEvent *ev ) {
   struct fy_event *fye;
   const char *str;

   memset( ev, 0, sizeof( *ev ) );
   ev->type = YAML_NO_EVENT;

   if( !parser->fyp ) {
      return 0;
   }

   fye = fy_parser_parse( parser->fyp );

   if( !fye ) {
      return 0;
   }

   ev->fye = fye;
   ev->fyp = parser->fyp;
   ev->type = fye->type;

   switch( fye->type ) {
   case FYET_SCALAR:
      if( fye->scalar.anchor ) {
         ev->anchor = fy_token_get_text0( fye->scalar.anchor );
      }
      if( fye->scalar.tag ) {
         ev->tag = fy_token_get_text0( fye->scalar.tag );
      }
      if( fye->scalar.value ) {
         str = fy_token_get_text0( fye->scalar.value );
         ev->value = str;
         ev->value_length = str ? strlen( str ) : 0;
      }
      break;
   case FYET_MAPPING_START:
      if( fye->mapping_start.anchor ) {
         ev->anchor = fy_token_get_text0( fye->mapping_start.anchor );
      }
      if( fye->mapping_start.tag ) {
         ev->tag = fy_token_get_text0( fye->mapping_start.tag );
      }
      break;
   case FYET_SEQUENCE_START:
      if( fye->sequence_start.anchor ) {
         ev->anchor = fy_token_get_text0( fye->sequence_start.anchor );
      }
      if( fye->sequence_start.tag ) {
         ev->tag = fy_token_get_text0( fye->sequence_start.tag );
      }
      break;
   case FYET_ALIAS:
      if( fye->alias.anchor ) {
         ev->anchor = fy_token_get_text0( fye->alias.anchor );
      }
      break;
   default:
      break;
   }

   return 1;
}

/* Event delete: use the fyp pointer stored in the event to free it. */
static inline void astYamlEventDelete( AstYamlEvent *ev ) {
   if( ev->fye && ev->fyp ) {
      fy_parser_event_free( ev->fyp, ev->fye );
      ev->fye = NULL;
      ev->fyp = NULL;
   }
}

/* Emitter lifecycle. */

static inline int astYamlEmitterInitialize( AstYamlEmitter *emitter ) {
   memset( emitter, 0, sizeof( *emitter ) );
   return 1;
}

static inline void astYamlEmitterDelete( AstYamlEmitter *emitter ) {
   if( emitter->fye ) {
      fy_emitter_destroy( emitter->fye );
      emitter->fye = NULL;
   }
   if( emitter->line_buf_len > 0 && emitter->write_cb ) {
      emitter->write_cb( emitter->write_data, emitter->line_buf,
                         emitter->line_buf_len );
   }
   free( emitter->line_buf );
   emitter->line_buf = NULL;
   emitter->line_buf_len = 0;
   emitter->line_buf_cap = 0;
}

static inline void astYamlEmitterSetOutput( AstYamlEmitter *emitter,
                                            AstYamlWriteCb cb,
                                            void *data ) {
   emitter->write_cb = cb;
   emitter->write_data = data;
}

static inline void astYamlEmitterSetIndent( AstYamlEmitter *emitter,
                                            int indent ) {
   emitter->indent = indent;
}

/* _astYamlEmitterCreate: create fy_emitter lazily on first use.
   NOTE: FYECF_INDENT is a macro that packs the indent level into the
   cfg flags. */
static inline int _astYamlEmitterCreate( AstYamlEmitter *emitter ) {
   struct fy_emitter_cfg cfg;

   if( emitter->fye ) {
      return 1;
   }

   memset( &cfg, 0, sizeof( cfg ) );
   cfg.output = _astFyamlWriteAdapter;
   cfg.userdata = emitter;
   cfg.flags = FYECF_DEFAULT | FYECF_INDENT( emitter->indent );
   emitter->fye = fy_emitter_create( &cfg );
   if( emitter->fye ) {
      struct fy_diag *diag = fy_emitter_get_diag( emitter->fye );
      if( diag ) {
         fy_diag_set_collect_errors( diag, true );
         fy_diag_unref( diag );
      }
   }
   return emitter->fye != NULL;
}

/* Translate AstYamlStyle to libfyaml scalar/node style enums. */
static inline enum fy_scalar_style _astFyamlToScalarStyle( AstYamlStyle s ) {
   switch( s ) {
   case AST_YAML_STYLE_PLAIN:
      return FYSS_PLAIN;
   case AST_YAML_STYLE_SINGLE_QUOTED:
      return FYSS_SINGLE_QUOTED;
   case AST_YAML_STYLE_DOUBLE_QUOTED:
      return FYSS_DOUBLE_QUOTED;
   default:
      return FYSS_ANY;
   }
}

static inline enum fy_node_style _astFyamlToNodeStyle( AstYamlStyle s ) {
   switch( s ) {
   case AST_YAML_STYLE_BLOCK:
      return FYNS_BLOCK;
   case AST_YAML_STYLE_FLOW:
      return FYNS_FLOW;
   default:
      return FYNS_ANY;
   }
}

/* Emit functions.
   NOTE: the variadic arguments to fy_emit_event_create differ by event
   type.  Calling conventions are:
     FYET_DOCUMENT_START: (int implicit, version*, tags*)
     FYET_SCALAR: (scalar_style, value, len, anchor, tag)
     FYET_MAPPING_START: (node_style, anchor, tag)
     FYET_SEQUENCE_START: (node_style, anchor, tag)
   anchor and tag are plain const char * (no separate length arg). */

static inline int astYamlEmitStreamStart( AstYamlEmitter *emitter ) {
   struct fy_event *fye;

   if( !_astYamlEmitterCreate( emitter ) ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_STREAM_START );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitStreamEnd( AstYamlEmitter *emitter ) {
   struct fy_event *fye;

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_STREAM_END );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

/* NOTE: tag directives (tag_handle/tag_prefix) are not forwarded to
   libfyaml; libfyaml uses fully-expanded tags internally and we always
   pass fully-qualified tag URIs from yamlchan.c anyway. */
static inline int astYamlEmitDocumentStart( AstYamlEmitter *emitter,
                                            int with_version,
                                            int with_tags,
                                            const char *tag_handle,
                                            const char *tag_prefix ) {
   struct fy_event *fye;
   (void) with_version; /* not used */
   (void) with_tags; /* not used */
   (void) tag_handle; /* not used */
   (void) tag_prefix; /* not used */

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_DOCUMENT_START,
                               0, NULL, NULL );
   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitDocumentEnd( AstYamlEmitter *emitter,
                                          int implicit ) {
   struct fy_event *fye;

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_DOCUMENT_END,
                               (bool) implicit );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitScalar( AstYamlEmitter *emitter,
                                     const char *anchor,
                                     const char *tag, const char *value,
                                     size_t len,
                                     int plain_implicit, int quoted_implicit,
                                     AstYamlStyle style ) {
   struct fy_event *fye;
   (void) plain_implicit; /* not used */
   (void) quoted_implicit; /* not used */

   if( !emitter->fye ) {
      return 0;
   }

   /* libyaml omits the str tag automatically; libfyaml does not. */
   if( tag && strcmp( tag, YAML_STR_TAG ) == 0 )
      tag = NULL;

   fye = fy_emit_event_create( emitter->fye, FYET_SCALAR,
      _astFyamlToScalarStyle( style ),
      value, len,
      anchor, tag );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitMappingStart( AstYamlEmitter *emitter,
                                           const char *anchor,
                                           const char *tag, int implicit,
                                           AstYamlStyle style ) {
   struct fy_event *fye;
   (void) implicit; /* not used */

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_MAPPING_START,
      _astFyamlToNodeStyle( style ),
      anchor, tag );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitMappingEnd( AstYamlEmitter *emitter ) {
   struct fy_event *fye;

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_MAPPING_END );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitSequenceStart( AstYamlEmitter *emitter,
                                            const char *anchor,
                                            const char *tag, int implicit,
                                            AstYamlStyle style ) {
   struct fy_event *fye;
   (void) implicit; /* not used */

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_SEQUENCE_START,
      _astFyamlToNodeStyle( style ),
      anchor, tag );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

static inline int astYamlEmitSequenceEnd( AstYamlEmitter *emitter ) {
   struct fy_event *fye;

   if( !emitter->fye ) {
      return 0;
   }

   fye = fy_emit_event_create( emitter->fye, FYET_SEQUENCE_END );

   if( !fye ) {
      return 0;
   }

   return fy_emit_event( emitter->fye, fye ) == 0;
}

#endif /* defined( YAML ) / defined( FYAML ) */

#endif /* YAML_BACKEND_H_INCLUDED */
