/*  =========================================================================
    wap_client - Wallet Client API

    Copyright (c) the Contributors as noted in the AUTHORS file.
                                                                
    (insert license text here)                                  
    =========================================================================
*/

/*
@header
    Description of class for man page.
@discuss
    Detailed discussion of the class, if any.
@end
*/

#include "wap_classes.h"

//  Forward reference to method arguments structure
typedef struct _client_args_t client_args_t;

//  This structure defines the context for a client connection
typedef struct {
    //  These properties must always be present in the client_t
    //  and are set by the generated engine. The cmdpipe gets
    //  messages sent to the actor; the msgpipe may be used for
    //  faster asynchronous message flows.
    zsock_t *cmdpipe;           //  Command pipe to/from caller API
    zsock_t *msgpipe;           //  Message pipe to/from caller API
    zsock_t *dealer;            //  Socket to talk to server
    wap_proto_t *message;       //  Message to/from server
    client_args_t *args;        //  Arguments from methods
    
    //  Own properties
    int heartbeat_timer;        //  Timeout for heartbeats to server
    int retries;                //  How many heartbeats we've tried
} client_t;

//  Include the generated client engine
#include "wap_client_engine.inc"

//  Allocate properties and structures for a new client instance.
//  Return 0 if OK, -1 if failed

static int
client_initialize (client_t *self)
{
    //  We'll ping the server once per second
    self->heartbeat_timer = 1000;
    return 0;
}

//  Free properties and structures for a client instance

static void
client_terminate (client_t *self)
{
    //  Destroy properties here
}


//  ---------------------------------------------------------------------------
//  connect_to_server_endpoint
//

static void
connect_to_server_endpoint (client_t *self)
{
    if (zsock_connect (self->dealer, "%s", self->args->endpoint)) {
        engine_set_exception (self, bad_endpoint_event);
        zsys_warning ("could not connect to %s", self->args->endpoint);
    }
}


//  ---------------------------------------------------------------------------
//  set_client_identity
//

static void
set_client_identity (client_t *self)
{
    wap_proto_set_identity (self->message, self->args->identity);
}


//  ---------------------------------------------------------------------------
//  use_connect_timeout
//

static void
use_connect_timeout (client_t *self)
{
    engine_set_timeout (self, self->args->timeout);
}


//  ---------------------------------------------------------------------------
//  client_is_connected
//

static void
client_is_connected (client_t *self)
{
    self->retries = 0;
    engine_set_connected (self, true);
    engine_set_timeout (self, self->heartbeat_timer);
}


//  ---------------------------------------------------------------------------
//  check_if_connection_is_dead
//

static void
check_if_connection_is_dead (client_t *self)
{
    //  We send at most 3 heartbeats before expiring the server
    if (++self->retries >= 3) {
        engine_set_timeout (self, 0);
        engine_set_connected (self, false);
        engine_set_exception (self, exception_event);
    }
}


//  ---------------------------------------------------------------------------
//  prepare_blocks_command
//

static void
prepare_blocks_command (client_t *self)
{
    wap_proto_set_block_ids (self->message, &self->args->block_ids);
    wap_proto_set_start_height (self->message, self->args->start_height);
}


//  ---------------------------------------------------------------------------
//  signal_have_blocks_ok
//

static void
signal_have_blocks_ok (client_t *self)
{
    zmsg_t *msg = wap_proto_get_msg_data (self->message);
    assert(msg != 0);
    zsock_send (self->cmdpipe, "s488p", "BLOCKS OK", wap_proto_status(self->message),
                wap_proto_start_height (self->message),
                wap_proto_curr_height (self->message),
                msg);
}


// ---------------------------------------------------------------------------
// prepare_start_mining_command
//

static void
prepare_start_mining_command (client_t *self)
{
    wap_proto_set_address (self->message, &self->args->address);
    wap_proto_set_thread_count (self->message, self->args->thread_count);
}

//  ---------------------------------------------------------------------------
//  prepare_send_raw_transaction_command
//

static void
prepare_send_raw_transaction_command (client_t *self)
{
    wap_proto_set_tx_as_hex (self->message, &self->args->tx_as_hex);
}


//  ---------------------------------------------------------------------------
//  signal_have_send_raw_transaction_ok
//

static void
signal_have_send_raw_transaction_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SEND RAW TRANSACTION OK",
                wap_proto_status (self->message));
}


//  ---------------------------------------------------------------------------
//  prepare_get_tx_command
//

static void
prepare_get_tx_command (client_t *self)
{
    wap_proto_set_tx_id (self->message, &self->args->tx_id);
    wap_proto_set_as_json (self->message, self->args->as_json);
}


//  ---------------------------------------------------------------------------
//  signal_have_get_tx_ok
//

static void
signal_have_get_tx_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p1", "GET TX OK",
                wap_proto_status (self->message),
                wap_proto_get_tx_data (self->message),
                wap_proto_in_pool (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_get_height_ok
//

static void
signal_have_get_height_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s48", "GET HEIGHT OK",
                wap_proto_status (self->message),
                wap_proto_height (self->message));
}


//  ---------------------------------------------------------------------------
//  signal_have_save_ok
//

static void
signal_have_save_bc_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SAVE BC OK", wap_proto_status(self->message));
}


//  ---------------------------------------------------------------------------
//  signal_have_start_mining_ok
//

static void
signal_have_start_mining_ok (client_t *self)
{
    zsock_send(self->cmdpipe, "s4", "START MINING OK", wap_proto_status(self->message));
}


//  ---------------------------------------------------------------------------
//  signal_have_stop_mining_ok
//

static void
signal_have_stop_mining_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "STOP MINING OK", wap_proto_status(self->message));
}


//  ---------------------------------------------------------------------------
//  signal_success
//

static void
signal_success (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SUCCESS", wap_proto_status(self->message));
}


//  ---------------------------------------------------------------------------
//  signal_bad_endpoint
//

static void
signal_bad_endpoint (client_t *self)
{
    zsock_send (self->cmdpipe, "sis", "FAILURE", -1, "Bad server endpoint");
}


//  ---------------------------------------------------------------------------
//  signal_failure
//

static void
signal_failure (client_t *self)
{
    zsock_send (self->cmdpipe, "sis", "FAILURE", -1, wap_proto_reason (self->message));
}


//  ---------------------------------------------------------------------------
//  check_status_code
//

static void
check_status_code (client_t *self)
{
    if (wap_proto_status (self->message) == WAP_PROTO_COMMAND_INVALID)
        engine_set_next_event (self, command_invalid_event);
    else
        engine_set_next_event (self, other_event);
}


//  ---------------------------------------------------------------------------
//  signal_unhandled_error
//

static void
signal_unhandled_error (client_t *self)
{
    zsys_error ("unhandled error code from server");
}


//  ---------------------------------------------------------------------------
//  signal_server_not_present
//

static void
signal_server_not_present (client_t *self)
{
    zsock_send (self->cmdpipe, "sis", "FAILURE", -1, "Server is not reachable");
}




//  ---------------------------------------------------------------------------
//  prepare_get_output_indexes_command
//

static void
prepare_get_output_indexes_command (client_t *self)
{
    wap_proto_set_tx_id (self->message, &self->args->tx_id);
}

//  ---------------------------------------------------------------------------
//  signal_have_output_indexes_ok
//

static void
signal_have_output_indexes_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p", "OUTPUT INDEXES OK",
        wap_proto_status (self->message),
        wap_proto_get_o_indexes (self->message));
}

void wap_client_test(bool verbose) {
}




//  ---------------------------------------------------------------------------
//  prepare_get_random_outs_command
//

static void
prepare_get_random_outs_command (client_t *self)
{
    wap_proto_set_amounts (self->message, &self->args->amounts);
    wap_proto_set_outs_count (self->message, self->args->outs_count);
}


//  ---------------------------------------------------------------------------
//  signal_have_random_outs_ok
//

static void
signal_have_random_outs_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p", "RANDOM OUTS OK",
        wap_proto_status (self->message),
        wap_proto_get_random_outputs (self->message));
}



//  ---------------------------------------------------------------------------
//  signal_have_get_info_ok
//

static void
signal_have_get_info_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s488888888881", "GET INFO OK",
        wap_proto_status (self->message),
        wap_proto_height (self->message),
        wap_proto_target_height (self->message),
        wap_proto_difficulty (self->message),
        wap_proto_tx_count (self->message),
        wap_proto_tx_pool_size (self->message),
        wap_proto_alt_blocks_count (self->message),
        wap_proto_outgoing_connections_count (self->message),
        wap_proto_incoming_connections_count (self->message),
        wap_proto_white_peerlist_size (self->message),
        wap_proto_grey_peerlist_size (self->message),
        wap_proto_testnet (self->message));
}


//  ---------------------------------------------------------------------------
//  signal_have_get_peer_list_ok
//

static void
signal_have_get_peer_list_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4pp", "GET PEER LIST OK",
        wap_proto_status (self->message),
        wap_proto_get_white_list (self->message),
        wap_proto_get_gray_list (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_get_mining_ok
//

static void
signal_have_get_mining_status_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4188p", "GET MINING STATUS OK",
        wap_proto_status (self->message),
        wap_proto_active (self->message),
        wap_proto_speed (self->message),
        wap_proto_thread_count (self->message),
        wap_proto_get_address (self->message));
}

//  ---------------------------------------------------------------------------
//  prepare_set_hash_log_rate_command
//

static void
prepare_set_log_hash_rate_command (client_t *self)
{
		wap_proto_set_visible (self->message, self->args->visible);
}

//  ---------------------------------------------------------------------------
//  signal_have_set_log_hash_rate_ok
//

static void
signal_have_set_log_hash_rate_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SET LOG HASH RATE OK",
        wap_proto_status (self->message));
}

//  ---------------------------------------------------------------------------
//  prepare_set_log_level_command
//

static void
prepare_set_log_level_command (client_t *self)
{
		wap_proto_set_level (self->message, self->args->level);
}

//  ---------------------------------------------------------------------------
//  signal_have_set_log_level_ok
//

static void
signal_have_set_log_level_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SET LOG LEVEL OK",
        wap_proto_status (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_start_save_graph_ok
//

static void
signal_have_start_save_graph_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "START SAVE GRAPH OK",
        wap_proto_status (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_stop_save_graph_ok
//

static void
signal_have_stop_save_graph_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "STOP SAVE GRAPH OK",
        wap_proto_status (self->message));
}

//  ---------------------------------------------------------------------------
//  prepare_get_block_hash_command
//

static void
prepare_get_block_hash_command (client_t *self)
{
		wap_proto_set_height (self->message, self->args->height);
}

//  ---------------------------------------------------------------------------
//  signal_have_get_block_hash_ok
//

static void
signal_have_get_block_hash_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p", "GET BLOCK HASH OK",
        wap_proto_status (self->message), wap_proto_get_hash (self->message));
}

//  ---------------------------------------------------------------------------
//  prepare_get_block_template_command
//

static void
prepare_get_block_template_command (client_t *self)
{
		wap_proto_set_reserve_size (self->message, self->args->reserve_size);
		wap_proto_set_address (self->message, &self->args->address);
}

//  ---------------------------------------------------------------------------
//  signal_have_get_block_template_ok
//

static void
signal_have_get_block_template_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4888pp", "GET BLOCK TEMPLATE OK",
        wap_proto_status (self->message), 
        wap_proto_reserved_offset (self->message),
        wap_proto_height (self->message),
        wap_proto_difficulty (self->message),
        wap_proto_get_prev_hash (self->message),
        wap_proto_get_block_template_blob (self->message));
}



//  ---------------------------------------------------------------------------
//  signal_have_get_hard_fork_info_ok
//

static void
signal_have_get_hard_fork_info_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s41144414", "GET HARD FORK INFO OK",
        wap_proto_status (self->message), 
        wap_proto_hfversion (self->message), 
        wap_proto_enabled (self->message), 
        wap_proto_window (self->message), 
        wap_proto_votes (self->message), 
        wap_proto_threshold (self->message), 
        wap_proto_voting (self->message), 
        wap_proto_hfstate (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_get_connections_list_ok
//

static void
signal_have_get_connections_list_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p", "GET CONNECTIONS LIST OK",
        wap_proto_status (self->message), 
        wap_proto_get_connections (self->message));
}

//  ---------------------------------------------------------------------------
//  prepare_stop_daemon_command
//

static void
prepare_stop_daemon_command (client_t *self)
{
    wap_proto_set_fast(self->message, self->args->fast);
}


//
//  ---------------------------------------------------------------------------
//  signal_have_stop_daemon_ok
//

static void
signal_have_stop_daemon_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "STOP DAEMON OK",(uint64_t)0);
}


//  ---------------------------------------------------------------------------
//  prepare_get_block_by_height_command
//

static void
prepare_get_block_by_height_command (client_t *self)
{
    wap_proto_set_height (self->message, self->args->height);
    wap_proto_set_header_only (self->message, self->args->header_only);
    wap_proto_set_as_json (self->message, self->args->as_json);
}


//  ---------------------------------------------------------------------------
//  signal_have_get_block_by_height_ok
//

static void
signal_have_get_block_by_height_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p118p4188p88", "GET BLOCK BY HEIGHT OK",
        wap_proto_status (self->message), 
        wap_proto_get_block (self->message),
        wap_proto_major_version (self->message),
        wap_proto_minor_version (self->message),
        wap_proto_timestamp (self->message),
        wap_proto_get_prev_hash (self->message),
        wap_proto_nonce (self->message),
        wap_proto_orphan (self->message),
        wap_proto_height (self->message),
        wap_proto_depth (self->message),
        wap_proto_get_hash (self->message),
        wap_proto_difficulty (self->message),
        wap_proto_reward (self->message));
}


//  ---------------------------------------------------------------------------
//  prepare_get_block_by_hash_command
//

static void
prepare_get_block_by_hash_command (client_t *self)
{
    wap_proto_set_hash (self->message, &self->args->hash);
    wap_proto_set_header_only (self->message, self->args->header_only);
    wap_proto_set_as_json (self->message, self->args->as_json);
}

//  ---------------------------------------------------------------------------
//  signal_have_get_block_by_hash_ok
//

static void
signal_have_get_block_by_hash_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p118p4188p88", "GET BLOCK BY HASH OK",
        wap_proto_status (self->message), 
        wap_proto_get_block (self->message),
        wap_proto_major_version (self->message),
        wap_proto_minor_version (self->message),
        wap_proto_timestamp (self->message),
        wap_proto_get_prev_hash (self->message),
        wap_proto_nonce (self->message),
        wap_proto_orphan (self->message),
        wap_proto_height (self->message),
        wap_proto_depth (self->message),
        wap_proto_get_hash (self->message),
        wap_proto_difficulty (self->message),
        wap_proto_reward (self->message));
}


//  ---------------------------------------------------------------------------
//  prepare_hard_fork_info_command
//

static void
prepare_hard_fork_info_command (client_t *self)
{
    wap_proto_set_hfversion(self->message, self->args->hfversion);
}


//  ---------------------------------------------------------------------------
//  prepare_get_key_image_status_command
//

static void
prepare_get_key_image_status_command (client_t *self)
{
    wap_proto_set_key_images(self->message, &self->args->key_images);
}


//  ---------------------------------------------------------------------------
//  signal_have_get_key_image_status_ok
//

static void
signal_have_get_key_image_status_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4p", "GET KEY IMAGE STATUS OK",
                wap_proto_status(self->message),
                wap_proto_get_spent (self->message));
}

//  ---------------------------------------------------------------------------
//  signal_have_get_tx_pool_ok
//

static void
signal_have_get_tx_pool_ok (client_t *self)
{
    zmsg_t *msg = wap_proto_get_msg_data (self->message);
    assert(msg != 0);
    zsock_send (self->cmdpipe, "s4p", "GET TX POOL OK",
                wap_proto_status(self->message),
                msg);
}

//  ---------------------------------------------------------------------------
//  prepare_set_out_peers_command
//

static void
prepare_set_out_peers_command (client_t *self)
{
    wap_proto_set_num_out_peers(self->message, self->args->num_out_peers);
}


//
//  ---------------------------------------------------------------------------
//  signal_have_set_out_peers_ok
//

static void
signal_have_set_out_peers_ok (client_t *self)
{
    zsock_send (self->cmdpipe, "s4", "SET OUT PEERS OK", wap_proto_status(self->message));
}


