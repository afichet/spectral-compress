# Debug tools for extracting framebuffers from dumps or files
add_subdirectory(exr-dump-extract)
add_subdirectory(jxl-dump-extract)
add_subdirectory(export-layers)
add_subdirectory(validate-image)

# Direct lossy compression without further processing to have a baseline
add_subdirectory(stupid-compress)
add_subdirectory(stupid-decompress)

# Image generation tool
add_subdirectory(generate-macbeth-img)

# Evaluation tools
add_subdirectory(quantization-from-image)
add_subdirectory(eval-quantization)

# Image conversion tools
add_subdirectory(cave-to-exr)
add_subdirectory(envi-to-exr)
add_subdirectory(exr-diff)
add_subdirectory(exr-png)

# REMOVE!
add_subdirectory(temp-diff)
