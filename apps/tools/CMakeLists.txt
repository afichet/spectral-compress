# Debug tools for extracting framebuffers from dumps or files
add_subdirectory(exr-dump-extract)
add_subdirectory(jxl-dump-extract)
add_subdirectory(export-layers)
add_subdirectory(validate-image)
add_subdirectory(sgeg-min-max)

# Direct lossy compression without further processing to have a baseline
add_subdirectory(simple-compress)
add_subdirectory(simple-decompress)

# Image generation tool
add_subdirectory(generate-macbeth-img)

# Evaluation tools
# add_subdirectory(quantization-from-image)
# add_subdirectory(eval-quantization)

# Image conversion tools
add_subdirectory(cave-exr)
add_subdirectory(envi-exr)
add_subdirectory(exr-png)
add_subdirectory(exr-pfm)

# Image visualization tools
add_subdirectory(exr-diff)
add_subdirectory(decompress-chrominance)
add_subdirectory(exr-peek-spectrum)

# Image manipulation tools
add_subdirectory(exr-strip-rgb)
add_subdirectory(exr-change-compression)
