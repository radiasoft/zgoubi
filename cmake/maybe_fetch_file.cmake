# Simple function to download assets and verify their checksum
function(maybe_fetch_file FILE URL SHA256 HAVE_FILE)
  get_filename_component(FILE_PATH ${FILE} DIRECTORY)
  if( (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${FILE}) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${FILE}) )
    message(STATUS "Copying ${FILE} from the source tree to the build tree.")
    file(COPY  ${CMAKE_CURRENT_SOURCE_DIR}/${FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${FILE_PATH})
  elseif(NOT EXISTS ${CMAKE_BINARY_DIR}/${FILE} )
    message( STATUS "Fetching ${FILE} from ${URL}")
    file(DOWNLOAD "${URL}" "${CMAKE_CURRENT_BINARY_DIR}/${FILE}"
      SHOW_PROGRESS
      TLS_VERIFY ON
      EXPECTED_HASH SHA256=${SHA256}
      )
    file(COPY ${CMAKE_CURRENT_BINARY_DIR}/${FILE} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/${FILE_PATH})
  endif()
  if( EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${FILE} )
    message( STATUS "Verifying SHA256 for ${FILE}")
    file(SHA256 ${CMAKE_CURRENT_BINARY_DIR}/${FILE} local_sha256)
    if( NOT ${local_sha256} STREQUAL ${SHA256} )
      message( FATAL_ERROR
	"Checksum for ${FILE} does not match expected value!
         ${CMAKE_CURRENT_BINARY_DIR}:${FILE} SHA256:
              ${local_sha256}
         expected SHA256:
              ${SHA256}")
      set(${HAVE_FILE} FALSE PARENT_SCOPE)
    else()
      message( STATUS "Checksum matches: ${FILE} SHA256=${local_sha256}")
      set(${HAVE_FILE} TRUE PARENT_SCOPE)
    endif()
  endif()
endfunction()
